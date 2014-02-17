{-# LANGUAGE GeneralizedNewtypeDeriving, TypeFamilies, FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses, BangPatterns #-}
-- | Common data types used everywhere.  This module is a collection of
-- very basic "bioinformatics" data types that are simple, but don't
-- make sense to define over and over.

module Bio.Base(
    Nucleotide(..),
    Qual(..), errorProb, fromErrorProb,
    ErrProb(..), toErrProb, fromErrProb, qualToErrProb,

    Word8,
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

import Bio.Util             ( phredplus, phredminus )
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

newtype Nucleotide = N { unN :: Word8 } deriving
    ( Eq, Ord, Ix, Storable
    , VG.Vector VU.Vector, VM.MVector VU.MVector, VU.Unbox )

instance Bounded Nucleotide where
    minBound = N  0
    maxBound = N 15

-- | Qualities are stored in deciban, also known as the Phred scale.  To
-- represent a value @p@, we store @-10 * log_10 p@.  Operations work
-- directly on the \"Phred\" value, as the name suggests.  The same goes
-- for the 'Ord' instance:  greater quality means higher \"Phred\"
-- score, meand lower error probability.
newtype Qual = Q { unQ :: Word8 } deriving
    ( Eq, Ord, Storable, Bounded, VG.Vector VU.Vector, VM.MVector VU.MVector, VU.Unbox )

instance Show Qual where
    showsPrec p (Q q) = (:) 'Q' . showsPrec p q

fromErrorProb :: (Floating a, RealFrac a) => a -> Qual
fromErrorProb a = Q $ round (-10 * log a / log 10)

errorProb :: Qual -> Double
errorProb (Q q) = 10 ** (-(fromIntegral q) / 10)

-- | A positive 'Double' value stored in log domain.  The scale is the
-- same \"Phred\" scale used for 'Qual' values, but here the semantics
-- derive from the stored 'Double' value.  In particular, the smaller
-- 'ErrProb' value would be the greater 'Qual'
newtype ErrProb = EP { unEP :: Double } deriving
    ( Eq, Storable, VG.Vector VU.Vector, VM.MVector VU.MVector, VU.Unbox )

instance Show ErrProb where
    showsPrec p (EP q) = (:) 'Q' . showsPrec p q

instance Ord ErrProb where
    EP a `compare` EP b = b `compare` a
    EP a   `min`   EP b = EP (a `max` b)
    EP a   `max`   EP b = EP (a `min` b)

    EP a <  EP b  =  b  < a
    EP a <= EP b  =  b <= a
    EP a >  EP b  =  b  > a
    EP a >= EP b  =  b >= a

instance Num ErrProb where
    fromInteger a = EP (-10 * log (fromInteger a) / log 10)
    EP a + EP b = EP (a `phredplus`  b)
    EP a - EP b = EP (a `phredminus` b)
    EP a * EP b = EP (a + b)
    negate    _ = error "no negative error probabilities"
    abs       x = x
    signum    _ = EP 0

instance Fractional ErrProb where
    fromRational a = EP (-10 * log (fromRational a) / log 10)
    EP a  /  EP b = EP (a - b)
    recip  (EP a) = EP (negate a)

toErrProb :: Double -> ErrProb
toErrProb p = EP (-10 * log p / log 10)

fromErrProb :: ErrProb -> Double
fromErrProb (EP q) = 10 ** (-q / 10)

qualToErrProb :: Qual -> ErrProb
qualToErrProb (Q q) = EP (fromIntegral q)

gap, nucA, nucC, nucG, nucT, nucN :: Nucleotide
gap  = N 0
nucA = N 1
nucC = N 2
nucG = N 4
nucT = N 8
nucN = N 15


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
-- have a negative sign and automatically make sense.  Position 0 could
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

{-# INLINE showNucleotide #-}
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
{-# INLINE revcompl #-}
revcompl :: [Nucleotide] -> [Nucleotide]
revcompl = reverse . map compl

-- | Complements a Nucleotide.
{-# INLINE compl #-}
compl :: Nucleotide -> Nucleotide
compl (N x) = N $ arr ! (x .&. 15)
  where
    arr :: UArray Word8 Word8
    !arr = listArray (0,15) [ 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 ]


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

