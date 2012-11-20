{-# LANGUAGE GeneralizedNewtypeDeriving, TypeFamilies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
-- | Common data types used everywhere.  This module is a collection of
-- very basic "bioinformatics" data types that are simple, but don't
-- make sense to define over and over.

module Bio.Base(
    Nucleotide(..),
    Sequence,
    nucA, nucC, nucG, nucT, nucN, gap,
    toNucleotide,
    showNucleotide,
    isGap,
    isBase,
    isProperBase,
    properBases,
    compl,
    revcompl,

    Seqid,
    unpackSeqid,
    packSeqid,
    shelve,

    Position(..),
    shiftPosition,
    p_is_reverse,

    Range(..),
    shiftRange,
    reverseRange,
    extendRange,
    insideRange,
    wrapRange,

    w2c,
    c2w,

    findAuxFile
) where

import Control.Monad        ( liftM )
import Data.Array.Unboxed
import Data.Bits
import Data.Char            ( isAlpha, isSpace, ord, toUpper )
import Data.Word            ( Word8 )
import Foreign.Storable     ( Storable(..) )
import System.Directory     ( doesFileExist )
import System.FilePath      ( (</>), isAbsolute, splitSearchPath )
import System.Environment   ( getEnvironment )

import qualified Data.ByteString.Char8 as S
import qualified Data.ByteString.Lazy as L

import qualified Data.Vector.Generic         as VG
import qualified Data.Vector.Generic.Mutable as VM
import qualified Data.Vector.Unboxed         as VU

import Data.ByteString.Internal ( c2w, w2c )

-- | A nucleotide base in an alignment.
-- Experience says we're dealing with Ns and gaps all the type, so
-- purity be damned, they are included as if they were real bases.
--
-- To allow @Nucleotide@s to be unpacked and incorparated into
-- containers, we choose to represent them the same way as the BAM file
-- format:  as a 4 bit wide field.  Gaps are encoded as 0 where they
-- make sense, N is 15.  

newtype Nucleotide = N { unN :: Word8 } deriving (Eq, Ord, Ix, Storable)

instance Bounded Nucleotide where
    minBound = N  0
    maxBound = N 15

gap, nucA, nucC, nucG, nucT, nucN :: Nucleotide
gap  = N 0
nucA = N 1
nucC = N 2
nucG = N 4
nucT = N 8
nucN = N 15

newtype instance VU.MVector s Nucleotide = MV_Nucleotide (VU.MVector s Word8)
newtype instance VU.Vector    Nucleotide = V_Nucleotide  (VU.Vector    Word8)

instance VM.MVector VU.MVector Nucleotide where
    {-# INLINE basicLength #-}
    basicLength (MV_Nucleotide v) = VM.basicLength v

    {-# INLINE basicUnsafeSlice #-}
    basicUnsafeSlice i l (MV_Nucleotide v) = MV_Nucleotide (VM.basicUnsafeSlice i l v)

    {-# INLINE basicOverlaps #-}
    basicOverlaps (MV_Nucleotide v) (MV_Nucleotide w) = VM.basicOverlaps v w 

    {-# INLINE basicUnsafeNew #-}
    basicUnsafeNew l = MV_Nucleotide `liftM` VM.basicUnsafeNew l
   
    {-# INLINE basicUnsafeRead #-}
    basicUnsafeRead (MV_Nucleotide v) i = N `liftM` VM.basicUnsafeRead v i
   
    {-# INLINE basicUnsafeWrite #-}
    basicUnsafeWrite (MV_Nucleotide v) i (N e) = VM.basicUnsafeWrite v i e


instance VG.Vector VU.Vector Nucleotide where
    {-# INLINE basicUnsafeFreeze #-}
    basicUnsafeFreeze (MV_Nucleotide v) = V_Nucleotide `liftM` VG.basicUnsafeFreeze v

    {-# INLINE basicUnsafeThaw #-}
    basicUnsafeThaw (V_Nucleotide v) = MV_Nucleotide `liftM` VG.basicUnsafeThaw v

    {-# INLINE basicLength #-}
    basicLength (V_Nucleotide v) = VG.basicLength v

    {-# INLINE basicUnsafeSlice #-}
    basicUnsafeSlice i l (V_Nucleotide v) = V_Nucleotide (VG.basicUnsafeSlice i l v) 

    {-# INLINE basicUnsafeIndexM #-}
    basicUnsafeIndexM (V_Nucleotide v) i = N `liftM` VG.basicUnsafeIndexM v i


instance VU.Unbox Nucleotide


-- | Sequence identifiers are ASCII strings
-- Since we tend to store them for a while, we use strict byte strings.
-- If you get a lazy bytestring from somewhere, use 'shelve' to convert
-- it for storage.  Use @unpackSeqid@ and @packSeqid@ to avoid the
-- import of @Data.ByteString@.
type Seqid = S.ByteString

-- | comparatively short Sequences 
type Sequence = VU.Vector Nucleotide

-- | Unpacks a @Seqid@ into a @String@
unpackSeqid :: Seqid -> String
unpackSeqid = S.unpack

-- | Packs a @String@ into a @Seqid@.  Only works for ASCII subset.
packSeqid :: String -> Seqid
packSeqid = S.pack

-- | Copies a lazy @L.ByteString@ into a strict @S.ByteString@.
-- A copy is forced, even if the lazy bytestring is a single chunk.
-- This makes sure bytestrings in long term storage don't hold onto
-- larger strings.
shelve :: L.ByteString -> S.ByteString
shelve s = case L.toChunks s of
    [ ] -> S.empty
    [c] -> S.copy c
    cs  -> S.concat cs

-- | Coordinates in a genome.
-- The position is zero-based, no questions about it.  Think of the
-- position as pointing to the crack between two bases: looking forward
-- you see the next base to the right, looking in the reverse direction
-- you see the complement of the first base to the left.  
--
-- To encode the strand, we (virtually) reverse-complement any sequence
-- and prepend it to the normal one.  That way, reversed coordinates
-- have a negative and automatically make sense.  Position 0 could
-- either be the beginning of the sequence or the end on the reverse
-- strand... that ambiguity shouldn't really matter.

data Position = Pos {
        p_seq   :: {-# UNPACK #-} !Seqid,   -- ^ sequence (e.g. some chromosome)
        p_start :: {-# UNPACK #-} !Int      -- ^ offset, zero-based
    } deriving (Show, Eq, Ord)

p_is_reverse :: Position -> Bool
p_is_reverse = (< 0) . p_start

-- | Ranges in genomes
-- We combine a position with a length.  In 'Range pos len', 'pos' is
-- always the start of a stretch of length 'len'.  Positions therefore
-- move in the opposite direction on the reverse strand.  To get the
-- same stretch on the reverse strand, shift r_pos by r_length, then
-- reverse direction (or call reverseRange).
data Range = Range {
        r_pos    :: {-# UNPACK #-} !Position, 
        r_length :: {-# UNPACK #-} !Int 
    } deriving (Show, Eq, Ord)


-- | Converts a character into a 'Nucleotide'.
-- The usual codes for A,C,G,T and U are understood, '-' and '.' become
-- gaps and everything else is an N.
toNucleotide :: Char -> Nucleotide
toNucleotide c = if inRange (bounds arr) (ord c) then N (arr ! ord c) else nucN
  where
    arr :: UArray Int Word8
    arr = listArray (0,127) (repeat (unN nucN)) //
          ( [ (ord x, unN n) | (x,n) <- pairs ] ++
            [ (ord (toUpper x), unN n) | (x,n) <- pairs ] )

    N a `plus` N b = N (a .|. b)

    pairs = [ ('a', nucA), ('c', nucC), ('g', nucG), ('t', nucT),
              ('u', nucT), ('-', gap),  ('.', gap),
              ('b', nucC `plus` nucG `plus` nucT), 
              ('d', nucA `plus` nucG `plus` nucT),
              ('h', nucA `plus` nucC `plus` nucT),
              ('v', nucA `plus` nucC `plus` nucG),
              ('k', nucG `plus` nucT),
              ('m', nucA `plus` nucC),
              ('s', nucC `plus` nucG), 
              ('w', nucA `plus` nucT),
              ('r', nucA `plus` nucG),
              ('y', nucC `plus` nucT) ]

-- | Tests if a 'Nucleotide' is a base.
-- Returns 'True' for everything but gaps.
isBase :: Nucleotide -> Bool
isBase (N n) = n /= 0

-- | Tests if a 'Nucleotide' is a proper base.
-- Returns 'True' for A,C,G,T only.
isProperBase :: Nucleotide -> Bool
isProperBase x = x == nucA || x == nucC || x == nucG || x == nucT

properBases :: [Nucleotide]
properBases = [ nucA, nucC, nucG, nucT ]

-- | Tests if a 'Nucleotide' is a gap.
-- Returns true only for the gap.
isGap :: Nucleotide -> Bool
isGap x = x == gap

showNucleotide :: Nucleotide -> Char
showNucleotide (N x) = S.index str $ fromIntegral $ x .&. 15
  where str = S.pack "-ACMGRSVTWYHKDBN"

instance Show Nucleotide where
    show     x = [ showNucleotide x ]
    showList l = (map showNucleotide l ++)

instance Read Nucleotide where
    readsPrec _ (c:cs) = [(toNucleotide c, cs)]
    readsPrec _ [    ] = []
    readList s = let (hd,tl) = span (\c -> isAlpha c || isSpace c || '-' == c) s
                 in [(map toNucleotide $ filter (not . isSpace) hd, tl)]
    
-- | Reverse-complements a stretch of Nucleotides
revcompl :: [Nucleotide] -> [Nucleotide]
revcompl = reverse . map compl

-- | Complements a Nucleotide.
compl :: Nucleotide -> Nucleotide
compl (N x) = N $ arr ! (x .&. 15)
  where
    arr :: UArray Word8 Word8
    arr = listArray (0,15) [ 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 ]


-- | Moves a @Position@.  The position is moved forward according to the
-- strand, negative indexes move backward accordingly.
shiftPosition :: Int -> Position -> Position
shiftPosition a p = p { p_start = p_start p + a }

-- | Moves a @Range@.  This is just @shiftPosition@ lifted.
shiftRange :: Int -> Range -> Range
shiftRange a r = r { r_pos = shiftPosition a (r_pos r) }

-- | Reverses a 'Range' to give the same @Range@ on the opposite strand.
reverseRange :: Range -> Range
reverseRange (Range (Pos sq pos) len) = Range (Pos sq (-pos-len)) len

-- | Extends a range.  The length of the range is simply increased.
extendRange :: Int -> Range -> Range
extendRange a r = r { r_length = r_length r + a }

-- | Expands a subrange.
-- @(range1 `insideRange` range2)@ interprets @range1@ as a subrange of
-- @range2@ and computes its absolute coordinates.  The sequence name of
-- @range1@ is ignored.
insideRange :: Range -> Range -> Range
insideRange r1@(Range (Pos _ start1) length1) r2@(Range (Pos sq start2) length2)
    | start2 < 0         = reverseRange (insideRange r1 (reverseRange r2))
    | start1 <= length2  = Range (Pos sq (start2 + start1)) (min length1 (length2 - start1))
    | otherwise          = Range (Pos sq (start2 + length2)) 0


-- | Wraps a range to a region.  This simply normalizes the start
-- position to be in the interval '[0,n)', which only makes sense if the
-- @Range@ is to be mapped onto a circular genome.  This works on both
-- strands and the strand information is retained.
wrapRange :: Int -> Range -> Range
wrapRange n (Range (Pos sq s) l) = Range (Pos sq (s `mod` n)) l

-- | Finds a file by searching the environment variable BIOHAZARD like a
-- PATH.
findAuxFile :: FilePath -> IO FilePath
findAuxFile fn | isAbsolute fn = return fn
               | otherwise = loop . maybe ["."] splitSearchPath . lookup "BIOHAZARD" =<< getEnvironment
  where
    loop [    ] = return fn
    loop (p:ps) = do e <- doesFileExist $ p </> fn
                     if e then return $ p </> fn else loop ps

