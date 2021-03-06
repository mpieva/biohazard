{-# LANGUAGE TypeFamilies #-}

-- | Parsers and Printers for BAM and SAM.  We employ an @Iteratee@
-- interface, and we strive to support everything possible in BAM.  So
-- far, the implementation of the nucleotides is somewhat lacking:  we
-- do not have support for ambiguity codes, and the "=" symbol is not
-- understood.

module Bio.Bam.Rec (
    BamRaw,
    bamRaw,
    virt_offset,
    raw_data,

    BamRec(..),
    unpackBam,
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
    isAlternative,
    isExactIndex,
    type_mask,

    progressBam,
    Word32
) where

import Bio.Bam.Header
import Bio.Iteratee
import Bio.Prelude

import Control.Monad.Primitive      ( unsafePrimToPrim, unsafeInlineIO )
import Foreign.C.Types              ( CInt(..), CSize(..) )
import Foreign.Marshal.Alloc        ( alloca )

import qualified Data.ByteString                    as B
import qualified Data.ByteString.Char8              as S
import qualified Data.ByteString.Internal           as B
import qualified Data.ByteString.Unsafe             as B
import qualified Data.Vector.Generic                as V
import qualified Data.Vector.Generic.Mutable        as VM
import qualified Data.Vector.Storable               as VS
import qualified Data.Vector.Unboxed                as U


-- | Cigar line in BAM coding
-- Bam encodes an operation and a length into a single integer, we keep
-- those integers in an array.
data Cigar = !CigOp :* !Int deriving (Eq, Ord)
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

-- | Extracts the aligned length from a cigar line.
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

-- | A vector that packs two 'Nucleotides' into one byte, just like Bam does.
data Vector_Nucs_half a = Vector_Nucs_half !Int !Int !(ForeignPtr Word8)

-- | A mutable vector that packs two 'Nucleotides' into one byte, just like Bam does.
data MVector_Nucs_half s a = MVector_Nucs_half !Int !Int !(ForeignPtr Word8)

type instance V.Mutable Vector_Nucs_half = MVector_Nucs_half

instance V.Vector Vector_Nucs_half Nucleotides where
    {-# INLINE basicUnsafeFreeze #-}
    basicUnsafeFreeze (MVector_Nucs_half o l fp) = return $  Vector_Nucs_half o l fp
    {-# INLINE basicUnsafeThaw #-}
    basicUnsafeThaw    (Vector_Nucs_half o l fp) = return $ MVector_Nucs_half o l fp

    {-# INLINE basicLength #-}
    basicLength          (Vector_Nucs_half _ l  _) = l
    {-# INLINE basicUnsafeSlice #-}
    basicUnsafeSlice s l (Vector_Nucs_half o _ fp) = Vector_Nucs_half (o + s) l fp

    {-# INLINE basicUnsafeIndexM #-}
    basicUnsafeIndexM (Vector_Nucs_half o _ fp) i
        | even (o+i) = return . Ns $! (b `shiftR` 4) .&. 0xF
        | otherwise  = return . Ns $!  b             .&. 0xF
      where !b = unsafeInlineIO $ withForeignPtr fp $ \p -> peekByteOff p ((o+i) `shiftR` 1)

instance VM.MVector MVector_Nucs_half Nucleotides where
    {-# INLINE basicLength #-}
    basicLength          (MVector_Nucs_half _ l  _) = l
    {-# INLINE basicUnsafeSlice #-}
    basicUnsafeSlice s l (MVector_Nucs_half o _ fp) = MVector_Nucs_half (o + s) l fp

    {-# INLINE basicOverlaps #-}
    basicOverlaps (MVector_Nucs_half _ _ fp1) (MVector_Nucs_half _ _ fp2) = fp1 == fp2
    {-# INLINE basicUnsafeNew #-}
    basicUnsafeNew l = unsafePrimToPrim $ MVector_Nucs_half 0 l <$> mallocForeignPtrBytes ((l+1) `shiftR` 1)

    {-# INLINE basicInitialize #-}
    basicInitialize v@(MVector_Nucs_half o l fp)

        | even    o = do unsafePrimToPrim $ withForeignPtr fp $ \p ->
                            memset (plusPtr p (o `shiftR` 1)) 0 (fromIntegral $ l `shiftR` 1)
                         when (odd l) $ VM.basicUnsafeWrite v (l-1) (Ns 0)

        | otherwise = do when (odd o) $ VM.basicUnsafeWrite v 0 (Ns 0)
                         unsafePrimToPrim $ withForeignPtr fp $ \p ->
                            memset (plusPtr p ((o+1) `shiftR` 1)) 0 (fromIntegral $ (l-1) `shiftR` 1)
                         when (even l) $ VM.basicUnsafeWrite v (l-1) (Ns 0)


    {-# INLINE basicUnsafeRead #-}
    basicUnsafeRead (MVector_Nucs_half o _ fp) i
        | even (o+i) = liftM (Ns . (.&.) 0xF . (`shiftR` 4)) b
        | otherwise  = liftM (Ns . (.&.) 0xF               ) b
      where b = unsafePrimToPrim $ withForeignPtr fp $ \p -> peekByteOff p ((o+i) `shiftR` 1)

    {-# INLINE basicUnsafeWrite #-}
    basicUnsafeWrite (MVector_Nucs_half o _ fp) i (Ns x) =
        unsafePrimToPrim $ withForeignPtr fp $ \p -> do
            y <- peekByteOff p ((o+i) `shiftR` 1)
            let y' | even (o+i) = x `shiftL` 4 .|. y .&. 0x0F
                   | otherwise  = x            .|. y .&. 0xF0
            pokeByteOff p ((o+i) `shiftR` 1) y'

foreign import ccall unsafe "string.h memset" memset
    :: Ptr Word8 -> CInt -> CSize -> IO ()

instance Show (Vector_Nucs_half Nucleotides) where
    show = show . V.toList

-- | Bam record in its native encoding along with virtual address.
data BamRaw = BamRaw { virt_offset :: {-# UNPACK #-} !FileOffset
                     , raw_data    :: {-# UNPACK #-} !Bytes }

-- | Smart constructor.  Makes sure we got a at least a full record.
{-# INLINE bamRaw #-}
bamRaw :: FileOffset -> Bytes -> BamRaw
bamRaw o s = if good then BamRaw o s else error $ "broken BAM record " ++ show (S.length s, m) ++ show m
  where
    good | S.length s < 32 = False
         | otherwise       = S.length s >= sum m
    m = [ 32, l_rnm, l_seq, (l_seq+1) `div` 2, l_cig * 4 ]
    l_rnm = fromIntegral (B.unsafeIndex s  8) - 1
    l_cig = fromIntegral (B.unsafeIndex s 12)             .|. fromIntegral (B.unsafeIndex s 13) `shiftL`  8
    l_seq = fromIntegral (B.unsafeIndex s 16)             .|. fromIntegral (B.unsafeIndex s 17) `shiftL`  8 .|.
            fromIntegral (B.unsafeIndex s 18) `shiftL` 16 .|. fromIntegral (B.unsafeIndex s 19) `shiftL` 24

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

        getInt8 :: Num a => Int -> a
        getInt8  o = fromIntegral (B.unsafeIndex (raw_data br) o)

        getInt16 :: (Num a, Bits a) => Int -> a
        getInt16 o = fromIntegral (B.unsafeIndex (raw_data br) o) .|.
                     fromIntegral (B.unsafeIndex (raw_data br) $ o+1) `shiftL`  8

        getInt32 :: (Num a, Bits a) => Int -> a
        getInt32 o = fromIntegral (        B.unsafeIndex (raw_data br) $ o+0 )             .|.
                     fromIntegral (        B.unsafeIndex (raw_data br) $ o+1 ) `shiftL`  8 .|.
                     fromIntegral (        B.unsafeIndex (raw_data br) $ o+2 ) `shiftL` 16 .|.
                     fromIntegral (signed (B.unsafeIndex (raw_data br) $ o+3)) `shiftL` 24

        signed :: Word8 -> Int8
        signed = fromIntegral

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

data Ext = Int Int | Float Float | Text Bytes | Bin Bytes | Char Word8
         | IntArr (U.Vector Int) | FloatArr (U.Vector Float)
    deriving (Show, Eq, Ord)

{-# INLINE unpackExtensions #-}
unpackExtensions :: Bytes -> Extensions
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


isPaired, isProperlyPaired, isUnmapped, isMateUnmapped, isReversed,
    isMateReversed, isFirstMate, isSecondMate, isAuxillary, isFailsQC,
    isDuplicate, isTrimmed, isMerged, isAlternative, isExactIndex :: BamRec -> Bool

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

isTrimmed        = flip testBit 0 . extAsInt 0 "FF"
isMerged         = flip testBit 1 . extAsInt 0 "FF"
isAlternative    = flip testBit 2 . extAsInt 0 "FF"
isExactIndex     = flip testBit 3 . extAsInt 0 "FF"

type_mask :: Int
type_mask = flagFirstMate .|. flagSecondMate .|. flagPaired

extAsInt :: Int -> BamKey -> BamRec -> Int
extAsInt d nm br = case lookup nm (b_exts br) of Just (Int i) -> i ; _ -> d

extAsString :: BamKey -> BamRec -> Bytes
extAsString nm br = case lookup nm (b_exts br) of
    Just (Char c) -> B.singleton c
    Just (Text s) -> s
    _             -> B.empty

setQualFlag :: Char -> BamRec -> BamRec
setQualFlag c br = br { b_exts = updateE "ZQ" (Text s') $ b_exts br }
  where
    s  = extAsString "ZQ" br
    s' = if c `S.elem` s then s else c `S.cons` s

-- | A simple progress indicator that prints sequence id and position.
progressBam :: MonadIO m => String -> Refs -> Int -> (String -> IO ()) -> Enumeratee [BamRaw] [BamRaw] m a
progressBam = progressPos (\br -> case unpackBam br of b -> (b_rname b, b_pos b))

