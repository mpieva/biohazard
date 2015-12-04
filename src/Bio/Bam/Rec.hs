{-# LANGUAGE OverloadedStrings, PatternGuards, BangPatterns #-}
{-# LANGUAGE NoMonomorphismRestriction, FlexibleContexts, FlexibleInstances #-}
{-# LANGUAGE RecordWildCards, TypeFamilies, MultiParamTypeClasses #-}
{-# LANGUAGE TemplateHaskell #-}

-- | Parsers and Printers for BAM and SAM.  We employ an @Iteratee@
-- interface, and we strive to support everything possible in BAM.  So
-- far, the implementation of the nucleotides is somewhat lacking:  we
-- do not have support for ambiguity codes, and the "=" symbol is not
-- understood.

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
    encodeBamEntry,
    encodeSamEntry,
    encodeBamWith2,
    encodeBamRawWith2,
    pushBam,
    pushBamRaw,
    writeBamFile2,
    pipeBamOutput2,
    writeBamHandle2,

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

    Cigar(..),
    CigOp(..),
    alignedLength,

    Nucleotides(..), Vector_Nucs_half,
    Extensions, Ext(..),
    extAsInt, extAsString, setQualFlag,
    deleteE, insertE, updateE, adjustE,

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
import Bio.Iteratee.Builder

import Control.Monad
import Control.Monad.Primitive      ( unsafePrimToPrim )
import Control.Applicative
import Data.Attoparsec.ByteString   ( anyWord8 )
import Data.Binary.Builder          ( toLazyByteString )
import Data.Binary.Put
import Data.Bits                    ( Bits, testBit, shiftL, shiftR, (.&.), (.|.) )
import Data.ByteString              ( ByteString )
import Data.ByteString.Internal     ( accursedUnutterablePerformIO )
import Data.Char                    ( chr, ord, digitToInt )
import Data.Foldable                ( foldMap )
import Data.Int                     ( Int32, Int16, Int8 )
import Data.Ix
import Data.Monoid
import Data.String                  ( fromString )
import Data.Word                    ( Word32, Word16 )
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc        ( alloca )
import Foreign.Storable             ( peek, poke, peekByteOff, pokeByteOff, Storable(..) )
import System.IO
import System.IO.Unsafe             ( unsafeDupablePerformIO )

import qualified Control.Monad.Catch                as C
import qualified Data.Attoparsec.ByteString.Char8   as P
import qualified Data.ByteString                    as B
import qualified Data.ByteString.Char8              as S
import qualified Data.ByteString.Internal           as B
import qualified Data.ByteString.Lazy.Char8         as L
import qualified Data.ByteString.Unsafe             as B
import qualified Data.Foldable                      as F
import qualified Data.Iteratee                      as I
import qualified Data.Map.Strict                    as M
import qualified Data.Vector.Generic                as V
import qualified Data.Vector.Generic.Mutable        as VM
import qualified Data.Vector.Storable               as VS
import qualified Data.Vector.Unboxed                as U
import qualified Data.Sequence                      as Z


-- | Cigar line in BAM coding
-- Bam encodes an operation and a length into a single integer, we keep
-- those integers in an array.
data Cigar = !CigOp :* !Int
infix 9 :*

data CigOp = Mat | Ins | Del | Nop | SMa | HMa | Pad
    deriving ( Eq, Ord, Enum, Show, Bounded, Ix )

instance Show Cigar where
    showsPrec _ (op :* num) = shows num . (:) (S.index "MIDNSHP" (fromEnum op))

instance Storable Cigar where
    sizeOf    _ = 4
    alignment _ = 1

    peek p = do w0 <- peekByteOff p 0 :: IO Word8
                w1 <- peekByteOff p 1 :: IO Word8
                w2 <- peekByteOff p 2 :: IO Word8
                w3 <- peekByteOff p 3 :: IO Word8
                let w = fromIntegral w0 `shiftL`  0 .|.  fromIntegral w1 `shiftL`  8 .|.
                        fromIntegral w2 `shiftL` 16 .|.  fromIntegral w3 `shiftL` 24
                return $ toEnum (w .&. 0xf) :* shiftR w 4

    poke p (op :* num) = do pokeByteOff p 0 (fromIntegral $ shiftR w  0 :: Word8)
                            pokeByteOff p 1 (fromIntegral $ shiftR w  8 :: Word8)
                            pokeByteOff p 2 (fromIntegral $ shiftR w 16 :: Word8)
                            pokeByteOff p 3 (fromIntegral $ shiftR w 24 :: Word8)
        where
            w = fromEnum op .|. shiftL num 4

-- | extracts the aligned length from a cigar line
-- This gives the length of an alignment as measured on the reference,
-- which is different from the length on the query or the length of the
-- alignment.
{-# INLINE alignedLength #-}
alignedLength :: V.Vector v Cigar => v Cigar -> Int
alignedLength = V.foldl' (\a -> (a +) . l) 0
  where l (op :* n) = if op == Mat || op == Del || op == Nop then n else 0


-- | internal representation of a BAM record
data BamRec = BamRec {
        b_qname :: Seqid,
        b_flag  :: Int,
        b_rname :: Refseq,
        b_pos   :: Int,
        b_mapq  :: Qual,
        b_cigar :: VS.Vector Cigar,
        b_mrnm  :: Refseq,
        b_mpos  :: Int,
        b_isize :: Int,
        b_seq   :: Vector_Nucs_half Nucleotides,
        b_qual  :: VS.Vector Qual,
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
        b_cigar = VS.empty,
        b_mrnm  = invalidRefseq,
        b_mpos  = invalidPos,
        b_isize = 0,
        b_seq   = V.empty,
        b_qual  = VS.empty,
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
    wrap enee it' = enee (\hdr -> I.mapStream unpackBam (it' hdr)) >>= lift . run


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
      where !b = accursedUnutterablePerformIO $ withForeignPtr fp $
                        \p -> peekByteOff p ((o+i) `shiftR` 1)

instance VM.MVector MVector_Nucs_half Nucleotides where
    basicLength          (MVector_Nucs_half _ l  _) = l
    basicUnsafeSlice s l (MVector_Nucs_half o _ fp) = MVector_Nucs_half (o + s) l fp

    basicOverlaps (MVector_Nucs_half _ _ fp1) (MVector_Nucs_half _ _ fp2) = fp1 == fp2
    basicUnsafeNew l = unsafePrimToPrim $ MVector_Nucs_half 0 l <$> mallocForeignPtrBytes ((l+1) `shiftR` 1)

    basicUnsafeRead (MVector_Nucs_half o _ fp) i
        | even (o+i) = liftM (Ns . (.&.) 0xF . (`shiftR` 4)) b
        | otherwise  = liftM (Ns . (.&.) 0xF               ) b
      where b = unsafePrimToPrim $ withForeignPtr fp $ \p -> peekByteOff p ((o+i) `shiftR` 1)

    basicUnsafeWrite (MVector_Nucs_half o _ fp) i (Ns x) =
        unsafePrimToPrim $ withForeignPtr fp $ \p -> do
            y <- peekByteOff p ((o+i) `shiftR` 1)
            let y' | even (o+i) = x `shiftL` 4 .|. y .&. 0x0F
                   | otherwise  = x            .|. y .&. 0xF0
            pokeByteOff p ((o+i) `shiftR` 1) y'

instance Show (Vector_Nucs_half Nucleotides) where
    show = show . V.toList

{-# INLINE[1] unpackBam #-}
unpackBam :: BamRaw -> BamRec
unpackBam br = BamRec {
        b_rname =      Refseq $ getInt32  0,
        b_pos   =               getInt32  4,
        b_mapq  =           Q $ getInt8   9,
        b_flag  =               getInt16 14,
        b_mrnm  =      Refseq $ getInt32 20,
        b_mpos  =               getInt32 24,
        b_isize = fromIntegral (getInt32 28 :: Int32),

        b_qname = B.unsafeTake l_read_name $ B.unsafeDrop 32 $ raw_data br,
        b_cigar = VS.unsafeCast $ VS.unsafeFromForeignPtr fp (off0+off_c) (4*l_cigar),
        b_seq   = Vector_Nucs_half (2 * (off_s+off0)) l_seq fp,
        b_qual  = VS.unsafeCast $ VS.unsafeFromForeignPtr fp (off0+off_q) l_seq,

        b_exts  = unpackExtensions $ S.drop off_e $ raw_data br,
        b_virtual_offset = virt_offset br }
  where
        (fp, off0, _) = B.toForeignPtr $ raw_data br
        off_c =    33 + l_read_name
        off_s = off_c + 4 * l_cigar
        off_q = off_s + (l_seq + 1) `div` 2
        off_e = off_q +  l_seq

        l_read_name = getInt8   8 - 1
        l_seq       = getInt32 16
        l_cigar     = getInt16 12

        getInt8 :: (Num a, Bits a) => Int -> a
        getInt8  o = fromIntegral (B.unsafeIndex (raw_data br) o)

        getInt16 :: (Num a, Bits a) => Int -> a
        getInt16 o = fromIntegral (B.unsafeIndex (raw_data br) o) .|.
                     fromIntegral (B.unsafeIndex (raw_data br) $ o+1) `shiftL`  8

        getInt32 :: (Num a, Bits a) => Int -> a
        getInt32 o = fromIntegral (B.unsafeIndex (raw_data br) $ o+0)             .|.
                     fromIntegral (B.unsafeIndex (raw_data br) $ o+1) `shiftL`  8 .|.
                     fromIntegral (B.unsafeIndex (raw_data br) $ o+2) `shiftL` 16 .|.
                     fromIntegral (B.unsafeIndex (raw_data br) $ o+3) `shiftL` 24

-- | A collection of extension fields.  The key is actually only two @Char@s, but that proved impractical.
-- (Hmm... we could introduce a Key type that is a 16 bit int, then give
-- it an @instance IsString@... practical?)
type Extensions = [( BamKey, Ext )]

-- | Deletes all occurences of some extension field.
deleteE :: BamKey -> Extensions -> Extensions
deleteE k = filter ((/=) k . fst)

-- | Blindly inserts an extension field.  This can create duplicates
-- (and there is no telling how other tools react to that).
insertE :: BamKey -> Ext -> Extensions -> Extensions
insertE k v = (:) (k,v)

-- | Deletes all occurences of an extension field, then inserts it with
-- a new value.  This is safer than 'insertE', but also more expensive.
updateE :: BamKey -> Ext -> Extensions -> Extensions
updateE k v = insertE k v . deleteE k

-- | Adjusts a named extension by applying a function.
adjustE :: (Ext -> Ext) -> BamKey -> Extensions -> Extensions
adjustE _ _ [         ]             = []
adjustE f k ((k',v):es) | k  ==  k' = (k', f v) : es
                        | otherwise = (k',   v) : adjustE f k es

data Ext = Int Int | Float Float | Text ByteString | Bin ByteString | Char Word8
         | IntArr (U.Vector Int) | FloatArr (U.Vector Float)
    deriving (Show, Eq, Ord)

{-# INLINE unpackExtensions #-}
unpackExtensions :: ByteString -> Extensions
unpackExtensions = go
  where
    go s | S.length s < 4 = []
         | otherwise = let key = fromString [ S.index s 0, S.index s 1 ]
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

    getFloat s = unsafeDupablePerformIO $ alloca $ \buf ->
                 pokeByteOff buf 0 (getInt 'I' s :: Word32) >> peek buf


encodeBamEntry :: BamRec -> BamRaw
encodeBamEntry = bamRaw 0 . S.concat . L.toChunks . runPut . putEntry
  where
    putEntry b | V.length (b_seq b) == V.length (b_qual b) = do
                     putWord32le   $ unRefseq $ b_rname b
                     put_int_32    $ b_pos b
                     put_int_8     $ S.length (b_qname b) + 1
                     put_int_8     $ unQ (b_mapq b)
                     put_int_16    $ distinctBin (b_pos b) (alignedLength (b_cigar b))
                     put_int_16    $ VS.length $ b_cigar b
                     put_int_16    $ b_flag b
                     put_int_32    $ V.length $ b_seq b
                     putWord32le   $ unRefseq $ b_mrnm b
                     put_int_32    $ b_mpos b
                     put_int_32    $ b_isize b
                     putByteString $ b_qname b
                     putWord8 0
                     VS.mapM_ putWord8 (VS.unsafeCast $ b_cigar b :: VS.Vector Word8)
                     putSeq $ b_seq b
                     VS.mapM_ (putWord8 . unQ) $ b_qual b
                     forM_ (more_exts b) $ \(BamKey k,v) -> putWord16le k >> putValue v

    more_exts :: BamRec -> Extensions
    more_exts b = if xf /= 0 then x' else b_exts b
        where xf = b_flag b `shiftR` 16
              x' = insertE "FF" (Int xf) $ b_exts b

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
pipeRawSamOutput hdr = joinI $ mapStream unpackBam $ pipeSamOutput hdr

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
    fromFloat float = unsafeDupablePerformIO $ alloca $ \buf ->
                      pokeByteOff buf 0 float >> peek buf


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


extAsInt :: Int -> BamKey -> BamRec -> Int
extAsInt d nm br = case lookup nm (b_exts br) of Just (Int i) -> i ; _ -> d

extAsString :: BamKey -> BamRec -> ByteString
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
parseSamRec ref = mkBamRec
                  <$> word <*> num <*> (ref <$> word) <*> (subtract 1 <$> num)
                  <*> (Q <$> num) <*> (VS.fromList <$> cigar) <*> rnext <*> (subtract 1 <$> num)
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

    quals    = {-# SCC "parseSamRec/quals" #-} defaultQs <$ P.char '*' <* sep <|> bsToVec <$> word
        where
            defaultQs sq = VS.replicate (V.length sq) (Q 0xff)
            bsToVec qs _ = VS.fromList . map (Q . subtract 33) $ B.unpack qs

    cigar    = [] <$ P.char '*' <* sep <|>
               P.manyTill (flip (:*) <$> P.decimal <*> cigop) sep

    cigop    = P.choice $ zipWith (\c r -> r <$ P.char c) "MIDNSHP" [Mat,Ins,Del,Nop,SMa,HMa,Pad]
    exts     = ext `P.sepBy` sep
    ext      = (\a b v -> (fromString [a,b],v)) <$> P.anyChar <*> P.anyChar <*> (P.char ':' *> value)

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

    mkBamRec nm fl rn po mq cg rn' mp is sq qs' =
                BamRec nm fl rn po mq cg (rn' rn) mp is sq (qs' sq)

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
    sarr = conjoin ',' . map shows . U.toList

-- | Encodes BAM records straight into a dynamic buffer, the BGZF's it.
-- Should be fairly direct and perform well.
{-# INLINE encodeBamWith2 #-}
encodeBamWith2 :: MonadIO m => Int -> BamMeta -> Enumeratee [BamRec] B.ByteString m a
encodeBamWith2 lv meta = joinI . eneeBam . encodeBgzfWith lv
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

{-# INLINE encodeBamRawWith2 #-}
encodeBamRawWith2 :: MonadIO m => Int -> BamMeta -> Enumeratee [BamRaw] B.ByteString m a
encodeBamRawWith2 lv meta = joinI . eneeBam . encodeBgzfWith lv
  where
    eneeBam  = eneeCheckIfDone (\k -> mapChunks (foldMap pushBamRaw) . k $ Chunk pushHeader)

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

writeBamFile2 :: FilePath -> BamMeta -> Iteratee [BamRec] IO ()
writeBamFile2 fp meta =
    C.bracket (liftIO $ openBinaryFile fp WriteMode)
              (liftIO . hClose)
              (flip writeBamHandle2 meta)

pipeBamOutput2 :: BamMeta -> Iteratee [BamRec] IO ()
pipeBamOutput2 meta = encodeBamWith2 0 meta =$ mapChunksM_ (liftIO . S.hPut stdout)

writeBamHandle2 :: MonadIO m => Handle -> BamMeta -> Iteratee [BamRec] m ()
writeBamHandle2 hdl meta = encodeBamWith2 6 meta =$ mapChunksM_ (liftIO . S.hPut hdl)

{-# RULES
    "pushBam/unpackBam"     forall b . pushBam (unpackBam b) = pushBamRaw b
  #-}

{-# INLINE[1] pushBam #-}
pushBam :: BamRec -> Push
pushBam BamRec{..} = mconcat
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

