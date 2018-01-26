-- | Printers for BAM and SAM.  BAM is properly supported, SAM can be
-- piped to standard output.

module Bio.Bam.Writer (
    IsBamRec(..),
    encodeBamWith,

    packBam,
    writeBamFile,
    writeBamHandle,
    pipeBamOutput,
    pipeSamOutput
                      ) where

import Bio.Bam.Header
import Bio.Bam.Rec
import Bio.Iteratee
import Bio.Iteratee.Builder
import Bio.Prelude

import Data.ByteString.Builder      ( hPutBuilder, Builder, toLazyByteString )
import Data.ByteString.Internal     ( ByteString(..) )
import Data.ByteString.Lazy         ( foldrChunks )
import Foreign.Marshal.Alloc        ( alloca )
import System.IO                    ( openBinaryFile, IOMode(..) )

import qualified Control.Monad.Catch                as C
import qualified Data.ByteString                    as B
import qualified Data.ByteString.Char8              as S
import qualified Data.Vector.Generic                as V
import qualified Data.Vector.Storable               as VS
import qualified Data.Vector.Unboxed                as U
import qualified Data.Sequence                      as Z

-- | write in SAM format to stdout
-- This is useful for piping to other tools (say, AWK scripts) or for
-- debugging.  No convenience function to send SAM to a file exists,
-- because that's a stupid idea.
pipeSamOutput :: MonadIO m => BamMeta -> Iteratee [BamRec] m ()
pipeSamOutput meta = do liftIO . hPutBuilder stdout $ showBamMeta meta
                        mapStreamM_ $ \b -> liftIO . putStr $ encodeSamEntry (meta_refs meta) b "\n"

encodeSamEntry :: Refs -> BamRec -> String -> String
encodeSamEntry refs b = conjoin '\t' [
    unpck (b_qname b),
    shows (b_flag b .&. 0xffff),
    unpck (sq_name $ getRef refs $ b_rname b),
    shows (b_pos b + 1),
    shows (unQ $ b_mapq b),
    V.foldr ((.) . shows) id (b_cigar b),
    if isValidRefseq (b_mrnm b) && b_mrnm b == b_rname b
      then (:) '=' else unpck (sq_name $ getRef refs $ b_mrnm b),
    shows (b_mpos b + 1),
    shows (b_isize b),
    shows (V.toList $ b_seq b),
    (++)  (V.toList . V.map (chr . (+33) . fromIntegral . unQ) $ b_qual b) ] .
    foldr (\(k,v) f -> (:) '\t' . shows k . (:) ':' . extToSam v . f) id (b_exts b)
  where
    unpck = (++) . S.unpack
    conjoin c = foldr1 (\a f -> a . (:) c . f)

    extToSam (Int      i) = (:) 'i' . (:) ':' . shows i
    extToSam (Float    f) = (:) 'f' . (:) ':' . shows f
    extToSam (Text     t) = (:) 'Z' . (:) ':' . unpck t
    extToSam (Bin      x) = (:) 'H' . (:) ':' . tohex x
    extToSam (Char     c) = (:) 'A' . (:) ':' . (:) (w2c c)
    extToSam (IntArr   a) = (:) 'B' . (:) ':' . (:) 'i' . sarr a
    extToSam (FloatArr a) = (:) 'B' . (:) ':' . (:) 'f' . sarr a

    tohex = B.foldr (\c f -> w2d (c `shiftR` 4) . w2d (c .&. 0xf) . f) id
    w2d = (:) . S.index "0123456789ABCDEF" . fromIntegral
    sarr v = conjoin ',' . map shows $ U.toList v

class IsBamRec a where
    pushBam :: a -> BgzfTokens -> BgzfTokens

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
encodeBamWith :: (MonadIO m, IsBamRec r) => Int -> BamMeta -> Enumeratee [r] S.ByteString m ()
encodeBamWith lv meta = eneeBam ><> encodeBgzf lv
  where
    eneeBam  = eneeCheckIfDone (\k -> mapChunks (foldMap (Endo . pushBam)) . k $ Chunk pushHeader)

    pushHeader :: Endo BgzfTokens
    pushHeader = Endo $ TkString "BAM\1"
                      . TkSetMark                        -- the length byte
                      . pushBuilder (showBamMeta meta)
                      . TkEndRecord                      -- fills the length in
                      . TkWord32 (fromIntegral . Z.length $ meta_refs meta)
                      . appEndo (foldMap (Endo . pushRef) (meta_refs meta))

    pushRef :: BamSQ -> BgzfTokens -> BgzfTokens
    pushRef bs = TkWord32 (fromIntegral $ B.length (sq_name bs) + 1)
               . TkString (sq_name bs)
               . TkWord8 0
               . TkWord32 (fromIntegral $ sq_length bs)

    pushBuilder :: Builder -> BgzfTokens -> BgzfTokens
    pushBuilder b tk = foldrChunks TkString tk (toLazyByteString b)

{-# INLINE pushBamRaw #-}
pushBamRaw :: BamRaw -> BgzfTokens -> BgzfTokens
pushBamRaw = TkLnString . raw_data

-- | Writes BAM encoded stuff to a file.
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

-- | Writes BAM encoded stuff to a 'Handle'.
writeBamHandle :: (MonadIO m, IsBamRec r) => Handle -> BamMeta -> Iteratee [r] m ()
writeBamHandle hdl meta = encodeBamWith 6 meta =$ mapChunksM_ (liftIO . S.hPut hdl)

{-# RULES
    "pushBam/unpackBam"     forall b . pushBamRec (unpackBam b) = pushBamRaw b
  #-}

{-# INLINE[1] pushBamRec #-}
pushBamRec :: BamRec -> BgzfTokens -> BgzfTokens
pushBamRec BamRec{..} =
      TkSetMark
    . TkWord32 (unRefseq b_rname)
    . TkWord32 (fromIntegral b_pos)
    . TkWord8  (fromIntegral $ B.length b_qname + 1)
    . TkWord8  (unQ b_mapq)
    . TkWord16 (fromIntegral bin)
    . TkWord16 (fromIntegral $ VS.length b_cigar)
    . TkWord16 (fromIntegral b_flag)
    . TkWord32 (fromIntegral $ V.length b_seq)
    . TkWord32 (unRefseq b_mrnm)
    . TkWord32 (fromIntegral b_mpos)
    . TkWord32 (fromIntegral b_isize)
    . TkString b_qname
    . TkWord8 0
    . VS.foldr ((.) . TkWord8) id (VS.unsafeCast b_cigar :: VS.Vector Word8)
    . pushSeq b_seq
    . VS.foldr ((.) . TkWord8 . unQ) id b_qual
    . foldr ((.) . pushExt) id b_exts
    . TkEndRecord
  where
    bin = distinctBin b_pos (alignedLength b_cigar)

    pushSeq :: V.Vector vec Nucleotides => vec Nucleotides -> BgzfTokens -> BgzfTokens
    pushSeq v = case v V.!? 0 of
                    Nothing -> id
                    Just a  -> case v V.!? 1 of
                        Nothing -> TkWord8 (unNs a `shiftL` 4)
                        Just b  -> TkWord8 (unNs a `shiftL` 4 .|. unNs b) . pushSeq (V.drop 2 v)

    pushExt :: (BamKey, Ext) -> BgzfTokens -> BgzfTokens
    pushExt (BamKey k, e) = case e of
        Text  t -> common 'Z' . TkString t . TkWord8 0
        Bin   t -> common 'H' . TkString t . TkWord8 0
        Char  c -> common 'A' . TkWord8 c
        Float f -> common 'f' . TkWord32 (fromIntegral $ fromFloat f)

        Int i   -> case put_some_int (U.singleton i) of
                        (c,op) -> common c . op i

        IntArr  ia -> case put_some_int ia of
                        (c,op) -> common 'B' . TkWord8 (fromIntegral $ ord c)
                                  . TkWord32 (fromIntegral $ U.length ia-1)
                                  . U.foldr ((.) . op) id ia

        FloatArr fa -> common 'B' . TkWord8 (fromIntegral $ ord 'f')
                       . TkWord32 (fromIntegral $ U.length fa-1)
                       . U.foldr ((.) . TkWord32 . fromFloat) id fa
      where
        common :: Char -> BgzfTokens -> BgzfTokens
        common z = TkWord16 k . TkWord8 (fromIntegral $ ord z)

        put_some_int :: U.Vector Int -> (Char, Int -> BgzfTokens -> BgzfTokens)
        put_some_int is
            | U.all (between        0    0xff) is = ('C', TkWord8  . fromIntegral)
            | U.all (between   (-0x80)   0x7f) is = ('c', TkWord8  . fromIntegral)
            | U.all (between        0  0xffff) is = ('S', TkWord16 . fromIntegral)
            | U.all (between (-0x8000) 0x7fff) is = ('s', TkWord16 . fromIntegral)
            | U.all                      (> 0) is = ('I', TkWord32 . fromIntegral)
            | otherwise                           = ('i', TkWord32 . fromIntegral)

        between :: Int -> Int -> Int -> Bool
        between l r x = l <= x && x <= r

        fromFloat :: Float -> Word32
        fromFloat float = unsafeDupablePerformIO $ alloca $ \buf ->
                          pokeByteOff buf 0 float >> peek buf

packBam :: BamRec -> IO BamRaw
packBam br = do bb <- newBuffer 1000
                (bb', TkEnd) <- store_loop bb (pushBamRec br TkEnd)
                return . bamRaw 0 $ PS (buffer bb') 4 (used bb' - 4)
  where
    store_loop bb tk = do (bb',tk') <- fillBuffer bb tk
                          case tk' of TkEnd -> return (bb',tk')
                                      _     -> do bb'' <- expandBuffer 0 bb'
                                                  store_loop bb'' tk'

