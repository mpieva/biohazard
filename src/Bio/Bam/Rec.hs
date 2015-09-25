{-# LANGUAGE OverloadedStrings, PatternGuards, BangPatterns #-}
{-# LANGUAGE NoMonomorphismRestriction, FlexibleContexts, FlexibleInstances #-}
{-# LANGUAGE RecordWildCards, TypeFamilies, MultiParamTypeClasses #-}

-- TODO:
-- - Automatic creation of some kind of index.  If possible, this should
--   be the standard index for sorted BAM and/or the newer CSI format.
--   Optionally, a block index for slicing of large files, even unsorted
--   ones.  Maybe an index by name and an index for group-sorted files.
--   Sensible indices should be generated whenever a file is written.
-- - Same for statistics.  Something like "flagstats" could always be
--   written.  Actually, having @writeBamHandle@ return enhanced
--   flagstats as a result might be even better.
--
-- TONOTDO:
-- - Reader for gzipped/bzipped/bgzf'ed SAM.  Storing SAM is a bad idea,
--   so why would anyone ever want to compress, much less index it?

module Bio.Bam.Rec (
    Block,
    BamEnumeratee,
    isBamOrSam,

    unpackBam,
    decodeBamEntry,
    encodeBamEntry,
    encodeSamEntry,

    decodeSam,
    decodeSam',

    decodeAnyBamOrSam,
    decodeAnyBamOrSamFile,

    writeBamFile,
    writeBamHandle,
    pipeBamOutput,
    pipeRawSamOutput,
    pipeSamOutput,

    BamRec(..),
    nullBamRec,
    getMd,

    Nucleotides(..), Vector_Nucs_half,
    Extensions, Ext(..),
    extAsInt, extAsString, setQualFlag,
    deleteE, insertE, updateE,

    isPaired,
    isProperlyPaired,
    isUnmapped,
    isMateUnmapped,
    isReversed,
    isMateReversed,
    isFirstMate,
    isSecondMate,
    isAuxillary,
    isFailsQC,
    isDuplicate,
    isTrimmed,
    isMerged,
    type_mask,

    Word32
) where

import Bio.Base
import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Iteratee

import Control.Monad
import Control.Monad.Primitive      ( unsafePrimToPrim )
import Control.Applicative
import Data.Attoparsec.ByteString   ( anyWord8 )
import Data.Binary.Builder          ( toLazyByteString )
import Data.Binary.Get
import Data.Binary.Put
import Data.Bits                    ( testBit, shiftL, shiftR, (.&.), (.|.), complement )
import Data.ByteString              ( ByteString )
import Data.ByteString.Internal     ( inlinePerformIO )
import Data.Char                    ( ord, digitToInt )
import Data.Int                     ( Int32, Int16, Int8 )
import Data.Monoid                  ( mempty )
import Data.Vector.Unboxed          ( (!?) )
import Data.Word                    ( Word32, Word16 )
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc        ( alloca )
import Foreign.Ptr                  ( castPtr )
import Foreign.Storable             ( peek, poke, peekByteOff, pokeByteOff )
import System.IO
import System.IO.Unsafe             ( unsafePerformIO )

import qualified Data.Attoparsec.ByteString.Char8   as P
import qualified Data.ByteString                    as B
import qualified Data.ByteString.Char8              as S
import qualified Data.ByteString.Internal           as B
import qualified Data.ByteString.Lazy.Char8         as L
import qualified Data.Foldable                      as F
import qualified Data.Iteratee                      as I
import qualified Data.Map                           as M
import qualified Data.Vector.Generic                as V
import qualified Data.Vector.Generic.Mutable        as VM
import qualified Data.Vector.Unboxed                as U

-- ^ Parsers and Printers for BAM and SAM.  We employ an @Iteratee@
-- interface, and we strive to support everything possible in BAM.  So
-- far, the implementation of the nucleotides is somewhat lacking:  we
-- do not have support for ambiguity codes, and the "=" symbol is not
-- understood.

-- | internal representation of a BAM record
{- data BamRec = BamRec {
        b_qname :: {-# UNPACK #-} !Seqid,
        b_flag  :: {-# UNPACK #-} !Int,
        b_rname :: {-# UNPACK #-} !Refseq,
        b_pos   :: {-# UNPACK #-} !Int,
        b_mapq  :: {-# UNPACK #-} !Qual,
        b_cigar :: Cigar,
        b_mrnm  :: {-# UNPACK #-} !Refseq,
        b_mpos  :: {-# UNPACK #-} !Int,
        b_isize :: {-# UNPACK #-} !Int,
        b_seq   :: !(Vector_Nucs_half Nucleotides),
        b_qual  :: !ByteString,         -- ^ quality, may be empty
        b_exts  :: Extensions,
        b_virtual_offset :: {-# UNPACK #-} !FileOffset -- ^ virtual offset for indexing purposes
    } deriving Show -}

-- | internal representation of a BAM record
data BamRec = BamRec {
        b_qname :: Seqid,
        b_flag  :: Int,
        b_rname :: Refseq,
        b_pos   :: Int,
        b_mapq  :: Qual,
        b_cigar :: Cigar,
        b_mrnm  :: Refseq,
        b_mpos  :: Int,
        b_isize :: Int,
        b_seq   :: Vector_Nucs_half Nucleotides,
        b_qual  :: ByteString,         -- ^ quality, may be empty
        b_exts  :: Extensions,
        b_virtual_offset :: FileOffset -- ^ virtual offset for indexing purposes
    } deriving Show

nullBamRec :: BamRec
nullBamRec = BamRec {
        b_qname = S.empty,
        b_flag  = flagUnmapped,
        b_rname = invalidRefseq,
        b_pos   = invalidPos,
        b_mapq  = Q 0,
        b_cigar = Cigar [],
        b_mrnm  = invalidRefseq,
        b_mpos  = invalidPos,
        b_isize = 0,
        b_seq   = V.empty,
        b_qual  = S.empty,
        b_exts  = [],
        b_virtual_offset = 0
    }

getMd :: BamRec -> Maybe [MdOp]
getMd r = case lookup "MD" $ b_exts r of
    Just (Text mdfield) -> readMd mdfield
    Just (Char mdfield) -> readMd $ B.singleton mdfield
    _                   -> Nothing

type BamEnumeratee m b = Enumeratee' BamMeta ByteString [BamRec] m b

isBamOrSam :: MonadIO m => Iteratee ByteString m (BamEnumeratee m a)
isBamOrSam = maybe decodeSam wrap `liftM` isBam
  where
    wrap enee it' = enee (\hdr -> I.mapStream decodeBamEntry (it' hdr)) >>= lift . run


-- | Checks if a file contains BAM in any of the common forms,
-- then decompresses it appropriately.  If the stream doesn't contain
-- BAM at all, it is instead decoded as SAM.  Since SAM is next to
-- impossible to recognize reliably, we don't even try.  Any old junk is
-- decoded as SAM and will fail later.
decodeAnyBamOrSam :: MonadIO m => BamEnumeratee m a
decodeAnyBamOrSam it = isBamOrSam >>= \k -> k it

decodeAnyBamOrSamFile :: (MonadIO m, MonadMask m)
                      => FilePath -> (BamMeta -> Iteratee [BamRec] m a) -> m (Iteratee [BamRec] m a)
decodeAnyBamOrSamFile fn k = enumFileRandom defaultBufSize fn (decodeAnyBamOrSam k) >>= run

-- | A vector that packs two 'Nucleotides' into one byte, just like Bam does.
data Vector_Nucs_half a = Vector_Nucs_half !Int !Int !(ForeignPtr Word8)

-- | A mutable vector that packs two 'Nucleotides' into one byte, just like Bam does.
data MVector_Nucs_half s a = MVector_Nucs_half !Int !Int !(ForeignPtr Word8)

type instance V.Mutable Vector_Nucs_half = MVector_Nucs_half

instance V.Vector Vector_Nucs_half Nucleotides where
    basicUnsafeFreeze (MVector_Nucs_half o l fp) = return $  Vector_Nucs_half o l fp
    basicUnsafeThaw    (Vector_Nucs_half o l fp) = return $ MVector_Nucs_half o l fp

    basicLength          (Vector_Nucs_half _ l  _) = l
    basicUnsafeSlice s l (Vector_Nucs_half o _ fp) = Vector_Nucs_half (o + s) l fp

    basicUnsafeIndexM (Vector_Nucs_half o _ fp) i
        | even (o+i) = return . Ns $ (b `shiftR` 4) .&. 0xF
        | otherwise  = return . Ns $  b             .&. 0xF
      where b = inlinePerformIO $ withForeignPtr fp $ \p -> peekByteOff p ((o+i) `shiftR` 1)

instance VM.MVector MVector_Nucs_half Nucleotides where
    basicLength          (MVector_Nucs_half _ l  _) = l
    basicUnsafeSlice s l (MVector_Nucs_half o _ fp) = MVector_Nucs_half (o + s) l fp

    basicOverlaps (MVector_Nucs_half _ _ fp1) (MVector_Nucs_half _ _ fp2) = fp1 == fp2
    basicUnsafeNew l = unsafePrimToPrim $ MVector_Nucs_half 0 l <$> mallocForeignPtrBytes ((l+1) `shiftR` 1)

    basicUnsafeRead (MVector_Nucs_half o _ fp) i
        | even (o+i) = Ns . (.&.) 0xF . (`shiftR` 4) <$> b
        | otherwise  = Ns . (.&.) 0xF                <$> b
      where b = unsafePrimToPrim $ withForeignPtr fp $ \p -> peekByteOff p ((o+i) `shiftR` 1)

    basicUnsafeWrite (MVector_Nucs_half o _ fp) i (Ns x) =
        unsafePrimToPrim $ withForeignPtr fp $ \p -> do
            y <- peekByteOff p ((o+i) `shiftR` 1)
            let y' | even (o+i) = x `shiftL` 4 .|. y .&. 0x0F
                   | otherwise  = x            .|. y .&. 0xF0
            pokeByteOff p ((o+i) `shiftR` 1) y'

instance Show (Vector_Nucs_half Nucleotides) where
    show = show . V.toList

{-# INLINE unpackBam #-}
unpackBam :: BamRaw -> BamRec
unpackBam br = BamRec {
        b_qname = br_qname br,
        b_flag  = br_flag br,
        b_rname = br_rname br,
        b_pos   = br_pos br,
        b_mapq  = br_mapq br,
        b_cigar = Cigar [ br_cigar_at br i | i <- [0..br_n_cigar_op br-1] ],
        b_mrnm  = br_mrnm br,
        b_mpos  = br_mpos br,
        b_isize = br_isize br,
        b_seq   = Vector_Nucs_half (2 * (off_s+off0)) (br_l_seq br) fp,
        b_qual  = S.take (br_l_seq br) $ S.drop off_q $ raw_data br,
        b_exts  = unpackExtensions $ S.drop off_e $ raw_data br,
        b_virtual_offset = virt_offset br }
  where
        (!fp, !off0, _) = B.toForeignPtr (raw_data br)
        !off_s = sum [ 33, br_l_read_name br, 4 * br_n_cigar_op br ]
        !off_q = off_s + (br_l_seq br + 1) `div` 2
        !off_e = off_q +  br_l_seq br

{-# DEPRECATED decodeBamEntry "Keep BamRaw if you can, use unpackBam if you must." #-}
-- | Decodes a raw block into a @BamRec@.
decodeBamEntry :: BamRaw -> BamRec
decodeBamEntry br = case pushEndOfInput $ runGetIncremental go `pushChunk` raw_data br of
        Fail _ _ m -> error m
        Partial  _ -> error "incomplete BAM record"
        Done _ _ r -> r
  where
    go = do !rid       <- Refseq       <$> getWord32le
            !start     <- fromIntegral <$> getWord32le
            !namelen   <- fromIntegral <$> getWord8
            !mapq      <-            Q <$> getWord8
            !_bin      <-                  getWord16le
            !cigar_len <- fromIntegral <$> getWord16le
            !flag      <- fromIntegral <$> getWord16le
            !read_len  <- fromIntegral <$> getWord32le
            !mate_rid  <- Refseq       <$> getWord32le
            !mate_pos  <- fromIntegral <$> getWord32le
            !ins_size  <- fromIntegral <$> getWord32le
            !read_name <- S.init       <$> getByteString namelen
            !cigar     <- Cigar . map decodeCigar <$> replicateM cigar_len getWord32le
            !qry_seq   <- getByteString $ (read_len+1) `div` 2
            !qual <- (\qs -> if B.all (0xff ==) qs then B.empty else qs) <$> getByteString read_len
            !exts <- getExtensions []

            return $ BamRec read_name flag rid start mapq cigar
                            mate_rid mate_pos ins_size
                            (V.fromListN read_len $ expand qry_seq)
                            qual exts (virt_offset br)

    expand t = if S.null t then [] else let x = B.head t in Ns (x `shiftR` 4) : Ns (x .&. 0xf) : expand (B.tail t)

    decodeCigar c | cc <= fromEnum (maxBound :: CigOp) = (toEnum cc, cl)
                  | otherwise = error "unknown Cigar operation"
      where cc = fromIntegral c .&. 0xf; cl = fromIntegral c `shiftR` 4

-- | A collection of extension fields.  The key is actually only two @Char@s, but that proved impractical.
-- (Hmm... we could introduce a Key type that is a 16 bit int, then give
-- it an @instance IsString@... practical?)
type Extensions = [( String, Ext )]

-- | Deletes all occurences of some extension field.
deleteE :: String -> Extensions -> Extensions
deleteE k = filter ((/=) k . fst)

-- | Blindly inserts an extension field.  This can create duplicates
-- (and there is no telling how other tools react to that).
insertE :: String -> Ext -> Extensions -> Extensions
insertE k v = (:) (k,v)

-- | Deletes all occurences of an extension field, then inserts it with
-- a new value.  This is safer than 'insertE', but also more expensive.
updateE :: String -> Ext -> Extensions -> Extensions
updateE k v = insertE k v . deleteE k


data Ext = Int Int | Float Float | Text ByteString | Bin ByteString | Char Word8
         | IntArr (U.Vector Int) | FloatArr (U.Vector Float)
    deriving (Show, Eq, Ord)

{-# INLINE unpackExtensions #-}
unpackExtensions :: ByteString -> Extensions
unpackExtensions = go
  where
    go s | S.length s < 4 = []
         | otherwise = let key = [ S.index s 0, S.index s 1 ]
                       in case S.index s 2 of
                         'Z' -> case S.break (== '\0') (S.drop 3 s) of (l,r) -> (key, Text l) : go (S.drop 1 r)
                         'H' -> case S.break (== '\0') (S.drop 3 s) of (l,r) -> (key, Bin  l) : go (S.drop 1 r)
                         'A' -> (key, Char (B.index s 3)) : go (S.drop 4 s)
                         'B' -> let tp = S.index s 3
                                    n  = getInt 'I' (S.drop 4 s)
                                in case tp of
                                      'f' -> (key, FloatArr (U.fromListN (n+1) [ getFloat (S.drop i s) | i <- [8, 12 ..] ]))
                                             : go (S.drop (12+4*n) s)
                                      _   -> (key, IntArr (U.fromListN (n+1) [ getInt tp (S.drop i s) | i <- [8, 8 + size tp ..] ]))
                                             : go (S.drop (8 + size tp * (n+1)) s)
                         'f' -> (key, Float (getFloat (S.drop 3 s))) : go (S.drop 7 s)
                         tp  -> (key, Int  (getInt tp (S.drop 3 s))) : go (S.drop (3 + size tp) s)

    size 'C' = 1
    size 'c' = 1
    size 'S' = 2
    size 's' = 2
    size 'I' = 4
    size 'i' = 4
    size 'f' = 4
    size  _  = 0

    getInt 'C' s | S.length s >= 1 = fromIntegral (fromIntegral (B.index s 0) :: Word8)
    getInt 'c' s | S.length s >= 1 = fromIntegral (fromIntegral (B.index s 0) ::  Int8)
    getInt 'S' s | S.length s >= 2 = fromIntegral                               (i :: Word16)
        where i = fromIntegral (B.index s 0) .|. fromIntegral (B.index s 1) `shiftL` 8
    getInt 's' s | S.length s >= 2 = fromIntegral                               (i ::  Int16)
        where i = fromIntegral (B.index s 0) .|. fromIntegral (B.index s 1) `shiftL` 8
    getInt 'I' s | S.length s >= 4 = fromIntegral                               (i :: Word32)
        where i = fromIntegral (B.index s 0)             .|. fromIntegral (B.index s 1) `shiftL`  8 .|.
                  fromIntegral (B.index s 2) `shiftL` 16 .|. fromIntegral (B.index s 3) `shiftL` 24
    getInt 'i' s | S.length s >= 4 = fromIntegral                               (i ::  Int32)
        where i = fromIntegral (B.index s 0)             .|. fromIntegral (B.index s 1) `shiftL`  8 .|.
                  fromIntegral (B.index s 2) `shiftL` 16 .|. fromIntegral (B.index s 3) `shiftL` 24
    getInt _ _ = 0

    getFloat s = unsafePerformIO $ alloca $ \buf ->
                 poke (castPtr buf) (getInt 'I' s :: Word32) >> peek buf


getExtensions :: Extensions -> Get Extensions
getExtensions m = getExt <|> return m
  where
    getExt :: Get Extensions
    getExt = do
            key <- (\a b -> [w2c a, w2c b]) <$> getWord8 <*> getWord8
            typ <- getWord8
            let cont v = getExtensions $ insertE key v m
            case w2c typ of
                    'Z' -> cont . Text =<< getByteStringNul
                    'H' -> cont . Bin  =<< getByteStringNul
                    'A' -> cont . Char =<< getWord8
                    'f' -> cont . Float . to_float =<< getWord32le
                    'B' -> do tp <- getWord8
                              n <- fromIntegral <$> getWord32le
                              case w2c tp of
                                 'f' -> cont . FloatArr . U.fromListN (n+1) . map to_float =<< replicateM (n+1) getWord32le
                                 x | Just get <- M.lookup x get_some_int -> cont . IntArr . U.fromListN (n+1) =<< replicateM (n+1) get
                                   | otherwise                           -> fail $ "array type code " ++ show x ++ " not recognized"
                    x | Just get <- M.lookup x get_some_int -> cont . Int =<< get
                      | otherwise                           -> fail $ "type code " ++ show x ++ " not recognized"

    to_float :: Word32 -> Float
    to_float word = unsafePerformIO $ alloca $ \buf ->
                    poke (castPtr buf) word >> peek buf

    get_some_int :: M.Map Char (Get Int)
    get_some_int = M.fromList $ zip "CcSsIi" [
                        fromIntegral                                     <$> getWord8,
                        fromIntegral . (fromIntegral ::  Word8 ->  Int8) <$> getWord8,
                        fromIntegral                                     <$> getWord16le,
                        fromIntegral . (fromIntegral :: Word16 -> Int16) <$> getWord16le,
                        fromIntegral                                     <$> getWord32le,
                        fromIntegral . (fromIntegral :: Word32 -> Int32) <$> getWord32le ]



getByteStringNul :: Get ByteString
getByteStringNul = S.init <$> (lookAhead (get_len 1) >>= getByteString)
  where
    get_len l = getWord8 >>= \w -> if w == 0 then return l else get_len $! l+1




encodeBamEntry :: BamRec -> BamRaw
encodeBamEntry = bamRaw 0 . S.concat . L.toChunks . runPut . putEntry
  where
    putEntry  b = do putWord32le   $ unRefseq $ b_rname b
                     put_int_32    $ b_pos b
                     put_int_8     $ S.length (b_qname b) + 1
                     put_int_8     $ unQ (b_mapq b)
                     put_int_16    $ distinctBin (b_pos b) (cigarToAlnLen (b_cigar b))
                     put_int_16    $ length $ unCigar $ b_cigar b
                     put_int_16    $ b_flag b
                     put_int_32    $ V.length $ b_seq b
                     putWord32le   $ unRefseq $ b_mrnm b
                     put_int_32    $ b_mpos b
                     put_int_32    $ b_isize b
                     putByteString $ b_qname b
                     putWord8 0
                     mapM_ (put_int_32 . encodeCigar) $ unCigar $ b_cigar b
                     putSeq $ b_seq b
                     putByteString $ if not (S.null (b_qual b)) then b_qual b
                                     else B.replicate (V.length $ b_seq b) 0xff
                     forM_ (more_exts b) $ \(k,v) ->
                        case k of [c,d] -> putChr c >> putChr d >> putValue v
                                  _     -> error $ "invalid field key " ++ show k

    more_exts :: BamRec -> Extensions
    more_exts b = if xf /= 0 then x' else b_exts b
        where xf = b_flag b `shiftR` 16
              x' = insertE "FF" (Int xf) $ b_exts b

    encodeCigar :: (CigOp,Int) -> Int
    encodeCigar (op,l) = fromEnum op .|. l `shiftL` 4

    putSeq v = case v V.!? 0 of
                 Nothing -> return ()
                 Just a  -> case v V.!? 1 of
                    Nothing -> putWord8 (unNs a `shiftL` 4)
                    Just b  -> do putWord8 (unNs a `shiftL` 4 .|. unNs b)
                                  putSeq (V.drop 2 v)

-- | writes BAM encoded stuff to a @Handle@
-- We generate BAM with dynamic blocks, then stream them out to the file.
--
-- XXX This could write indexes on the side---a simple block index
-- for MapReduce style slicing, a standard BAM index or a name index
-- would be possible.
writeBamHandle :: Handle -> BamMeta -> Iteratee [BamRec] IO ()
writeBamHandle hdl meta = I.mapStream encodeBamEntry =$ writeRawBamHandle hdl meta

-- | writes BAM encoded stuff to a file
-- XXX This should(!) write indexes on the side---a simple block index
-- for MapReduce style slicing, a standard BAM index or a name index
-- would be possible.  When writing to a file, this makes even more
-- sense than when writing to a @Handle@.
writeBamFile :: FilePath -> BamMeta -> Iteratee [BamRec] IO ()
writeBamFile fp meta = I.mapStream encodeBamEntry =$ writeRawBamFile fp meta

-- | write BAM encoded stuff to stdout
-- This send uncompressed BAM to stdout.  Useful for piping to other
-- tools.
pipeBamOutput :: BamMeta -> Iteratee [BamRec] IO ()
pipeBamOutput meta = I.mapStream encodeBamEntry =$ pipeRawBamOutput meta

-- | write in SAM format to stdout
-- This is useful for piping to other tools (say, AWK scripts) or for
-- debugging.  No convenience function to send SAM to a file exists,
-- because that's a stupid idea.
pipeSamOutput :: MonadIO m => BamMeta -> Iteratee [BamRec] m ()
pipeSamOutput meta = do liftIO . L.putStr . toLazyByteString $ showBamMeta meta
                        mapStreamM_ $ \b -> liftIO . putStr $ encodeSamEntry (meta_refs meta) b "\n"

pipeRawSamOutput :: MonadIO m => BamMeta -> Iteratee [BamRaw] m ()
pipeRawSamOutput hdr = joinI $ mapStream decodeBamEntry $ pipeSamOutput hdr

put_int_32, put_int_16, put_int_8 :: Integral a => a -> Put
put_int_32 = putWord32le . fromIntegral
put_int_16 = putWord16le . fromIntegral
put_int_8  = putWord8 . fromIntegral

putChr :: Char -> Put
putChr = putWord8 . fromIntegral . ord

putValue :: Ext -> Put
putValue v = case v of
    Text t      -> putChr 'Z' >> putByteString t >> putWord8 0
    Bin b       -> putChr 'H' >> putByteString b >> putWord8 0
    Char c      -> putChr 'A' >> putWord8 c
    Float f     -> putChr 'f' >> put_int_32 (fromFloat f)
    Int i       -> case put_some_int (U.singleton i) of
                        (c,op) -> putChr c >> op i
    IntArr   ia -> case put_some_int ia of
                        (c,op) -> putChr 'B' >> putChr c >> put_int_32 (U.length ia-1)
                                  >> mapM_ op (U.toList ia)
    FloatArr fa -> putChr 'B' >> putChr 'f' >> put_int_32 (U.length fa-1)
                   >> mapM_ (put_int_32 . fromFloat) (U.toList fa)
  where
    put_some_int :: U.Vector Int -> (Char, Int -> Put)
    put_some_int is
        | U.all (between        0    0xff) is = ('C', put_int_8)
        | U.all (between   (-0x80)   0x7f) is = ('c', put_int_8)
        | U.all (between        0  0xffff) is = ('S', put_int_16)
        | U.all (between (-0x8000) 0x7fff) is = ('s', put_int_16)
        | U.all                      (> 0) is = ('I', put_int_32)
        | otherwise                           = ('i', put_int_32)

    between :: Int -> Int -> Int -> Bool
    between l r x = l <= x && x <= r

    fromFloat :: Float -> Int32
    fromFloat float = unsafePerformIO $ alloca $ \buf ->
                      poke (castPtr buf) float >> peek buf


isPaired, isProperlyPaired, isUnmapped, isMateUnmapped, isReversed,
    isMateReversed, isFirstMate, isSecondMate, isAuxillary, isFailsQC,
    isDuplicate, isTrimmed, isMerged :: BamRec -> Bool

isPaired         = flip testBit  0 . b_flag
isProperlyPaired = flip testBit  1 . b_flag
isUnmapped       = flip testBit  2 . b_flag
isMateUnmapped   = flip testBit  3 . b_flag
isReversed       = flip testBit  4 . b_flag
isMateReversed   = flip testBit  5 . b_flag
isFirstMate      = flip testBit  6 . b_flag
isSecondMate     = flip testBit  7 . b_flag
isAuxillary      = flip testBit  8 . b_flag
isFailsQC        = flip testBit  9 . b_flag
isDuplicate      = flip testBit 10 . b_flag
isTrimmed        = flip testBit 16 . b_flag
isMerged         = flip testBit 17 . b_flag


type_mask :: Int
type_mask = flagFirstMate .|. flagSecondMate .|. flagPaired


extAsInt :: Int -> String -> BamRec -> Int
extAsInt d nm br = case lookup nm (b_exts br) of Just (Int i) -> i ; _ -> d

extAsString :: String -> BamRec -> ByteString
extAsString nm br = case lookup nm (b_exts br) of
    Just (Char c) -> B.singleton c
    Just (Text s) -> s
    _             -> B.empty


setQualFlag :: Char -> BamRec -> BamRec
setQualFlag c br = br { b_exts = updateE "ZQ" (Text s') $ b_exts br }
  where
    s  = extAsString "ZQ" br
    s' = if c `S.elem` s then s else c `S.cons` s


-- | Iteratee-style parser for SAM files, designed to be compatible with
-- the BAM parsers.  Parses plain uncompressed SAM, nothing else.  Since
-- it is supposed to work the same way as the BAM parser, it requires
-- the presense of the SQ header lines.  These are stripped from the
-- header text and turned into the symbol table.
decodeSam :: Monad m => (BamMeta -> Iteratee [BamRec] m a) -> Iteratee ByteString m (Iteratee [BamRec] m a)
decodeSam inner = joinI $ enumLinesBS $ do
    let pHeaderLine acc str = case P.parseOnly parseBamMetaLine str of Right f -> return $ f : acc
                                                                       Left e  -> fail $ e ++ ", " ++ show str
    meta <- liftM (foldr ($) mempty . reverse) (joinI $ I.breakE (not . S.isPrefixOf "@") $ I.foldM pHeaderLine [])
    decodeSamLoop (meta_refs meta) (inner meta)

decodeSamLoop :: Monad m => Refs -> Enumeratee [ByteString] [BamRec] m a
decodeSamLoop refs inner = I.convStream (liftI parse_record) inner
  where !refs' = M.fromList $ zip [ nm | BamSQ { sq_name = nm } <- F.toList refs ] [toEnum 0..]
        ref x = M.findWithDefault invalidRefseq x refs'

        parse_record (EOF x) = icont parse_record x
        parse_record (Chunk []) = liftI parse_record
        parse_record (Chunk (l:ls)) | "@" `S.isPrefixOf` l = parse_record (Chunk ls)
        parse_record (Chunk (l:ls)) = case P.parseOnly (parseSamRec ref) l of
            Right  r -> idone [r] (Chunk ls)
            Left err -> icont parse_record (Just $ iterStrExc $ err ++ ", " ++ show l)

-- | Parser for SAM that doesn't look for a header.  Has the advantage
-- that it doesn't stall on a pipe that never delivers data.  Has the
-- disadvantage that it never reads the header and therefore needs a
-- list of allowed RNAMEs.
decodeSam' :: Monad m => Refs -> Enumeratee ByteString [BamRec] m a
decodeSam' refs inner = joinI $ enumLinesBS $ decodeSamLoop refs inner

parseSamRec :: (ByteString -> Refseq) -> P.Parser BamRec
parseSamRec ref = (\nm fl rn po mq cg rn' -> BamRec nm fl rn po mq cg (rn' rn))
                  <$> word <*> num <*> (ref <$> word) <*> (subtract 1 <$> num)
                  <*> (Q <$> num) <*> (Cigar <$> cigar) <*> rnext <*> (subtract 1 <$> num)
                  <*> snum <*> sequ <*> quals <*> exts <*> pure 0
  where
    sep      = P.endOfInput <|> () <$ P.char '\t'
    word     = P.takeTill ((==) '\t') <* sep
    num      = P.decimal <* sep
    snum     = P.signed P.decimal <* sep

    rnext    = id <$ P.char '=' <* sep <|> const . ref <$> word
    sequ     = {-# SCC "parseSamRec/sequ" #-}
               (V.empty <$ P.char '*' <|>
               V.fromList . map toNucleotides . S.unpack <$> P.takeWhile is_nuc) <* sep

    quals    = {-# SCC "parseSamRec/quals" #-} B.empty <$ P.char '*' <* sep <|> B.map (subtract 33) <$> word

    cigar    = [] <$ P.char '*' <* sep <|>
               P.manyTill (flip (,) <$> P.decimal <*> cigop) sep

    cigop    = P.choice $ zipWith (\c r -> r <$ P.char c) "MIDNSHP" [Mat,Ins,Del,Nop,SMa,HMa,Pad]
    exts     = ext `P.sepBy` sep
    ext      = (\a b v -> ([a,b],v)) <$> P.anyChar <*> P.anyChar <*> (P.char ':' *> value)

    value    = P.char 'A' *> P.char ':' *> (Char <$>               anyWord8) <|>
               P.char 'i' *> P.char ':' *> (Int  <$>     P.signed P.decimal) <|>
               P.char 'Z' *> P.char ':' *> (Text <$> P.takeTill ((==) '\t')) <|>
               P.char 'H' *> P.char ':' *> (Bin  <$>               hexarray) <|>
               P.char 'f' *> P.char ':' *> (Float . realToFrac <$> P.double) <|>
               P.char 'B' *> P.char ':' *> (
                    P.satisfy (P.inClass "cCsSiI") *> (intArr   <$> many (P.char ',' *> P.signed P.decimal)) <|>
                    P.char 'f'                     *> (floatArr <$> many (P.char ',' *> P.double)))

    intArr   is = IntArr   $ U.fromList is
    floatArr fs = FloatArr $ U.fromList $ map realToFrac fs
    hexarray    = B.pack . repack . S.unpack <$> P.takeWhile (P.inClass "0-9A-Fa-f")
    repack (a:b:cs) = fromIntegral (digitToInt a * 16 + digitToInt b) : repack cs ; repack _ = []
    is_nuc = P.inClass "acgtswkmrybdhvnACGTSWKMRYBDHVN"

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
    unpck (B.map (+33) $ b_qual b) ] .
    foldr (\(k,v) f -> (:) '\t' . (++) k . (:) ':' . extToSam v . f) id (b_exts b)
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
    sarr = conjoin ',' . map shows . U.toList

