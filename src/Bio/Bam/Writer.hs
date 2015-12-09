{-# LANGUAGE RecordWildCards, OverloadedStrings, FlexibleContexts #-}
module Bio.Bam.Writer (
    IsBamRec(..),
    encodeBamWith,

    writeBamFile,
    writeBamHandle,
    pipeBamOutput,
    pipeSamOutput
                      ) where

import Bio.Base
import Bio.Bam.Header
import Bio.Bam.Rec
import Bio.Iteratee
import Bio.Iteratee.Builder

import Control.Applicative
import Data.ByteString.Builder      ( toLazyByteString )
import Data.Bits
import Data.Char                    ( ord, chr )
import Data.Monoid
import Foreign.Marshal.Alloc        ( alloca )
import Foreign.Storable             ( pokeByteOff, peek )
import System.IO
import System.IO.Unsafe             ( unsafeDupablePerformIO )

import qualified Control.Monad.Catch                as C
import qualified Data.ByteString                    as B
import qualified Data.ByteString.Char8              as S
import qualified Data.ByteString.Lazy               as L
import qualified Data.Vector.Generic                as V
import qualified Data.Vector.Storable               as VS
import qualified Data.Vector.Unboxed                as U
import qualified Data.Sequence                      as Z

-- ^ Printers for BAM.  We employ an @Iteratee@ interface, and we strive
-- to keep BAM records in their encoded form.  This is most compact and
-- often faster, since it saves the time for repeated decoding and
-- encoding, if that's not strictly needed.


-- | write in SAM format to stdout
-- This is useful for piping to other tools (say, AWK scripts) or for
-- debugging.  No convenience function to send SAM to a file exists,
-- because that's a stupid idea.
pipeSamOutput :: MonadIO m => BamMeta -> Iteratee [BamRec] m ()
pipeSamOutput meta = do liftIO . L.putStr . toLazyByteString $ showBamMeta meta
                        mapStreamM_ $ \b -> liftIO . putStr $ encodeSamEntry (meta_refs meta) b "\n"

encodeSamEntry :: Refs -> BamRec -> String -> String
encodeSamEntry refs b = conjoin '\t' [
    unpck (b_qname b),
    shows (b_flag b .&. 0xffff),
    unpck (sq_name $ getRef refs $ b_rname b),
    shows (b_pos b + 1),
    shows (b_mapq b),
    shows (b_cigar b),
    unpck (sq_name $ getRef refs $ b_mrnm b),
    shows (b_mpos b + 1),
    shows (b_isize b + 1),
    shows (V.toList $ b_seq b),
    (++)  (V.toList . V.map (chr . (+33) . fromIntegral . unQ) $ b_qual b) ] .
    foldr (\(k,v) f -> (:) '\t' . shows k . (:) ':' . extToSam v . f) id (b_exts b)
  where
    unpck = (++) . S.unpack
    conjoin c = foldr1 (\a f -> a . (:) c . f)

    extToSam (Int        i) = (:) 'i' . (:) ':' . shows i
    extToSam (Float      f) = (:) 'f' . (:) ':' . shows f
    extToSam (Text       t) = (:) 'Z' . (:) ':' . unpck t
    extToSam (Bin        x) = (:) 'H' . (:) ':' . tohex x
    extToSam (Char       c) = (:) 'A' . (:) ':' . (:) (w2c c)
    extToSam (IntArr   arr) = (:) 'B' . (:) ':' . (:) 'i' . sarr arr
    extToSam (FloatArr arr) = (:) 'B' . (:) ':' . (:) 'f' . sarr arr

    tohex = B.foldr (\c f -> w2d (c `shiftR` 4) . w2d (c .&. 0xf) . f) id
    w2d = (:) . S.index "0123456789ABCDEF" . fromIntegral
    sarr v = conjoin ',' . map shows $ U.toList v

class IsBamRec a where
    pushBam :: a -> Push

instance IsBamRec BamRaw where
    {-# INLINE pushBam #-}
    pushBam = pushBamRaw

instance IsBamRec BamRec where
    {-# INLINE pushBam #-}
    pushBam = pushBamRec

instance (IsBamRec a, IsBamRec b) => IsBamRec (Either a b) where
    {-# INLINE pushBam #-}
    pushBam = either pushBam pushBam

-- | Encodes BAM records straight into a dynamic buffer, the BGZF's it.
-- Should be fairly direct and perform well.
{-# INLINE encodeBamWith #-}
encodeBamWith :: (MonadIO m, IsBamRec r) => Int -> BamMeta -> Enumeratee [r] B.ByteString m a
encodeBamWith lv meta = joinI . eneeBam . encodeBgzfWith lv
  where
    eneeBam  = eneeCheckIfDone (\k -> mapChunks (foldMap pushBam) . k $ Chunk pushHeader)

    pushHeader = pushByteString "BAM\1"
              <> setMark                        -- the length byte
              <> pushBuilder (showBamMeta meta)
              <> endRecord                      -- fills the length in
              <> pushWord32 (fromIntegral . Z.length $ meta_refs meta)
              <> foldMap pushRef (meta_refs meta)

    pushRef bs = ensureBuffer     (fromIntegral $ B.length (sq_name bs) + 9)
              <> unsafePushWord32 (fromIntegral $ B.length (sq_name bs) + 1)
              <> unsafePushByteString (sq_name bs)
              <> unsafePushByte 0
              <> unsafePushWord32 (fromIntegral $ sq_length bs)

{-# INLINE pushBamRaw #-}
pushBamRaw :: BamRaw -> Push
pushBamRaw br = ensureBuffer (B.length (raw_data br) + 4)
             <> unsafePushWord32 (fromIntegral $ B.length (raw_data br))
             <> unsafePushByteString (raw_data br)

-- | writes BAM encoded stuff to a file
-- XXX This should(!) write indexes on the side---a simple block index
-- for MapReduce style slicing, a standard BAM index or a name index
-- would be possible.  When writing to a file, this makes even more
-- sense than when writing to a @Handle@.
writeBamFile :: IsBamRec r => FilePath -> BamMeta -> Iteratee [r] IO ()
writeBamFile fp meta =
    C.bracket (liftIO $ openBinaryFile fp WriteMode)
              (liftIO . hClose)
              (flip writeBamHandle meta)

-- | write BAM encoded stuff to stdout
-- This send uncompressed BAM to stdout.  Useful for piping to other
-- tools.
pipeBamOutput :: IsBamRec r => BamMeta -> Iteratee [r] IO ()
pipeBamOutput meta = encodeBamWith 0 meta =$ mapChunksM_ (liftIO . S.hPut stdout)

-- | writes BAM encoded stuff to a @Handle@
-- We generate BAM with dynamic blocks, then stream them out to the file.
--
-- XXX This could write indexes on the side---a simple block index
-- for MapReduce style slicing, a standard BAM index or a name index
-- would be possible.
writeBamHandle :: (MonadIO m, IsBamRec r) => Handle -> BamMeta -> Iteratee [r] m ()
writeBamHandle hdl meta = encodeBamWith 6 meta =$ mapChunksM_ (liftIO . S.hPut hdl)

{-# RULES
    "pushBam/unpackBam"     forall b . pushBamRec (unpackBam b) = pushBamRaw b
  #-}

{-# INLINE[1] pushBamRec #-}
pushBamRec :: BamRec -> Push
pushBamRec BamRec{..} = mconcat
    [ ensureBuffer minlength
    , unsafeSetMark
    , unsafePushWord32 $ unRefseq b_rname
    , unsafePushWord32 $ fromIntegral b_pos
    , unsafePushByte   $ fromIntegral $ B.length b_qname + 1
    , unsafePushByte   $ unQ b_mapq
    , unsafePushWord16 $ fromIntegral bin
    , unsafePushWord16 $ fromIntegral $ VS.length b_cigar
    , unsafePushWord16 $ fromIntegral b_flag
    , unsafePushWord32 $ fromIntegral $ V.length b_seq
    , unsafePushWord32 $ unRefseq b_mrnm
    , unsafePushWord32 $ fromIntegral b_mpos
    , unsafePushWord32 $ fromIntegral b_isize
    , unsafePushByteString b_qname
    , unsafePushByte 0
    , VS.foldr ((<>) . unsafePushByte) mempty (VS.unsafeCast b_cigar :: VS.Vector Word8)
    , pushSeq b_seq
    , VS.foldr ((<>) . unsafePushByte . unQ) mempty b_qual
    , foldMap pushExt b_exts
    , endRecord ]
  where
    bin = distinctBin b_pos (alignedLength b_cigar)
    minlength = 37 + B.length b_qname + 4 * V.length b_cigar + V.length b_qual + (V.length b_seq + 1) `shiftR` 1

    pushSeq :: V.Vector vec Nucleotides => vec Nucleotides -> Push
    pushSeq v = case v V.!? 0 of
                    Nothing -> mempty
                    Just a  -> case v V.!? 1 of
                        Nothing -> unsafePushByte (unNs a `shiftL` 4)
                        Just b  -> unsafePushByte (unNs a `shiftL` 4 .|. unNs b)
                                   <> pushSeq (V.drop 2 v)

    pushExt :: (BamKey, Ext) -> Push
    pushExt (BamKey k, e) = case e of
        Text t -> common (4 + B.length t) 'Z' $
                  unsafePushByteString t <> unsafePushByte 0

        Bin  t -> common (4 + B.length t) 'H' $
                  unsafePushByteString t <> unsafePushByte 0

        Char c -> common 4 'A' $ unsafePushByte c

        Float f -> common 7 'f' $ unsafePushWord32 (fromIntegral $ fromFloat f)

        Int i   -> case put_some_int (U.singleton i) of
                        (c,op) -> common 7 c (op i)

        IntArr  ia -> case put_some_int ia of
                        (c,op) -> common (4 * U.length ia) 'B' $ unsafePushByte (fromIntegral $ ord c)
                                  <> unsafePushWord32 (fromIntegral $ U.length ia-1)
                                  <> U.foldr ((<>) . op) mempty ia

        FloatArr fa -> common (4 * U.length fa) 'B' $ unsafePushByte (fromIntegral $ ord 'f')
                       <> unsafePushWord32 (fromIntegral $ U.length fa-1)
                       <> U.foldr ((<>) . unsafePushWord32 . fromFloat) mempty fa
      where
        common l z b = ensureBuffer l <> unsafePushWord16 k
                    <> unsafePushByte (fromIntegral $ ord z) <> b

        put_some_int :: U.Vector Int -> (Char, Int -> Push)
        put_some_int is
            | U.all (between        0    0xff) is = ('C', unsafePushByte . fromIntegral)
            | U.all (between   (-0x80)   0x7f) is = ('c', unsafePushByte . fromIntegral)
            | U.all (between        0  0xffff) is = ('S', unsafePushWord16 . fromIntegral)
            | U.all (between (-0x8000) 0x7fff) is = ('s', unsafePushWord16 . fromIntegral)
            | U.all                      (> 0) is = ('I', unsafePushWord32 . fromIntegral)
            | otherwise                           = ('i', unsafePushWord32 . fromIntegral)

        between :: Int -> Int -> Int -> Bool
        between l r x = l <= x && x <= r

        fromFloat :: Float -> Word32
        fromFloat float = unsafeDupablePerformIO $ alloca $ \buf ->
                          pokeByteOff buf 0 float >> peek buf

