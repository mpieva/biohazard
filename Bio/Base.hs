module Bio.Base(
    Nucleotide(..),
    toNucleotide,
    isGap,
    isBase,
    isProperBase,
    minBase,
    maxBase,
    compl,
    revcompl,

    Sense(..),

    Seqid,
    unpackSeqid,
    shelve,

    Position(..),
    shift_pos,

    Range(..),
    shift_range,
    reverse_range,
    extend,
    inside,
    wraprange
) where

import Data.Char            ( isAlpha, isSpace )
import Data.Ix              ( Ix )

import qualified Data.ByteString.Char8 as S
import qualified Data.ByteString.Lazy as L

-- Common data types used everywhere

-- | A base in an alignment.
-- Experience says we're dealing with Ns and gaps all the type, so
-- purity be damned, they are included as if they were real bases.
data Nucleotide = Gap | A | C | G | T | N deriving (Eq, Ord, Enum, Ix, Bounded)

-- | Sense of a strand.
-- Avoids the confusion inherent in using a simple bool.
data Sense = Forward | Reverse deriving (Show, Eq, Ord)

-- | Sequence identifiers are ACSII strings
-- Since we tend to store them for a while, we use strict byte strings.
-- If you get a lazu bytestring from somewhere, use 'shelve' to convert
-- it for storage.
type Seqid = S.ByteString

unpackSeqid :: Seqid -> String
unpackSeqid = S.unpack

-- | copies a lazy bytestring into a strict one
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
data Position = Pos {
        p_seq   :: {-# UNPACK #-} !Seqid,
        p_sense :: {-# UNPACK #-} !Sense,
        p_start :: {-# UNPACK #-} !Int
    } deriving (Show, Eq, Ord)

-- | Ranges in genomes
-- We combine a position with a length.  In 'Range pos len', 'pos' is
-- always the start of a stretch of length 'len'.  Positions therefore
-- move in the opposite direction on the reverse strand.  Shifting r_pos
-- by r_length, then reversing direction yields the same stretch on the
-- reverse strand.
data Range = Range {
        r_pos    :: {-# UNPACK #-} !Position, 
        r_length :: {-# UNPACK #-} !Int 
    } deriving (Show, Eq, Ord)


-- | converts a character into a 'Nucleotide'
-- The usual codes for A,C,G,T and U are understood, '-' and '.' become
-- gaps and everything else is an N.
toNucleotide :: Char -> Nucleotide
toNucleotide 'A' = A
toNucleotide 'C' = C
toNucleotide 'G' = G
toNucleotide 'T' = T
toNucleotide 'U' = T
toNucleotide 'a' = A
toNucleotide 'c' = C
toNucleotide 'g' = G
toNucleotide 't' = T
toNucleotide 'u' = T
toNucleotide '-' = Gap
toNucleotide '.' = Gap
toNucleotide  _  = N

-- | Tests if a 'Nucleotide' is a base.
-- Returns 'True' for everything but gaps.
isBase :: Nucleotide -> Bool
isBase Gap = False
isBase  _  = True

isProperBase :: Nucleotide -> Bool
isProperBase Gap = False
isProperBase  N  = False
isProperBase  _  = True

-- | Tests if a 'Nucleotide' is a gap.
-- Returns true only for the gap.
isGap :: Nucleotide -> Bool
isGap  Gap = True 
isGap   _  = False

instance Show Nucleotide where
    show Gap = "-"
    show  A  = "A"
    show  C  = "C"
    show  G  = "G"
    show  T  = "T"
    show  _  = "N"
    showList l = (concatMap show l ++)

instance Read Nucleotide where
    readsPrec _ (c:cs) = [(toNucleotide c, cs)]
    readsPrec _ [    ] = []
    readList s = let (hd,tl) = span (\c -> isAlpha c || isSpace c || '-' == c) s
                 in [(map toNucleotide $ filter (not . isSpace) hd, tl)]
    
-- | returns the smallest Nucleotide
-- This returns the smallest Nucleotide according to the Ord instance
-- that is not a gap.
minBase :: Nucleotide
minBase = A

-- | returns the largest Nucleotide
-- This returns the largest Nucleotide according to the Ord instance
-- that is not a gap.
maxBase :: Nucleotide
maxBase = N

-- | reverse-complements a stretch of Nucleotides
revcompl :: [Nucleotide] -> [Nucleotide]
revcompl = reverse . map compl

-- | complements a Nucleotide
compl :: Nucleotide -> Nucleotide
compl A = T
compl C = G
compl G = C
compl T = A
compl x = x


-- | moves a position
-- The position is moved forward according to the strand, negative
-- indexes move backward accordingly.
shift_pos :: Int -> Position -> Position
shift_pos a p = case p_sense p of Forward -> p { p_start = p_start p + a }
                                  Reverse -> p { p_start = p_start p - a }

shift_range :: Int -> Range -> Range
shift_range a r = r { r_pos = shift_pos a (r_pos r) }

reverse_range :: Range -> Range
reverse_range (Range (Pos sq Forward pos) len) = Range (Pos sq Reverse (pos+len)) len
reverse_range (Range (Pos sq Reverse pos) len) = Range (Pos sq Forward (pos-len)) len

-- | extends a range
-- The length of the range is simply increased.
extend :: Int -> Range -> Range
extend a r = r { r_length = r_length r + a }


-- | expands a subrange
-- (range1 `inside` range2) interprets range1 as a subrange of range2
-- and computes its absolute coordinates.  The sequence name range1 is
-- ignored.
inside :: Range -> Range -> Range
inside (Range (Pos  _ dir2 start2) length2) (Range (Pos sq dir1 start1) length1) =
    Range (Pos sq dir' pos') length2
  where
    dir' = if dir1 == dir2 then Forward else Reverse
    pos' = case dir1 of Forward -> start1 + start2
                        Reverse -> start1 - start2

-- | wraps a range to a ragion
-- This simply normalizes the start position to be in the interval [0,n).
wraprange :: Int -> Range -> Range
wraprange n (Range (Pos sq str s) l) = Range (Pos sq str (s `mod` n)) l

