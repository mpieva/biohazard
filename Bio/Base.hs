module Bio.Base(
    Nucleotide(..),
    toNucleotide,
    isGap,
    isBase,
    minBase,
    maxBase,
    compl,
    revcompl,

    Sense(..),

    Seqid,
    shelve,

    Position(..),
    shift,

    Range(..),
    extend
) where

import Data.Char            ( isAlpha, isSpace )

import qualified Data.ByteString      as S
import qualified Data.ByteString.Lazy as L

-- Common data types used everywhere

-- | A base in an alignment.
-- Experience says we're dealing with Ns and gaps all the type, so
-- purity be damned, they are included as if they were real bases.
data Nucleotide = Gap | A | C | G | T | N deriving (Eq, Ord, Enum)

-- | Sense of a strand.
-- Avoids the confusion inherent in using a simple bool.
data Sense = Forward | Reverse deriving (Show, Eq, Ord)

-- | Sequence identifiers are ACSII strings
-- Since we tend to store them for a while, we use strict byte strings.
-- If you get a lazu bytestring from somewhere, use 'shelve' to convert
-- it for storage.
type Seqid = S.ByteString

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
-- The position is zero-based, no questions about it.
data Position = Pos {
        p_seq :: !Seqid,
        p_sense :: !Sense,
        p_start :: !Int
    } deriving (Show, Eq, Ord)

-- | Ranges in genomes
-- We combine a position with a length.  In 'Range pos len', 'pos' is
-- always the start of a stretch of length 'len'.  Positions therefore
-- move in the opposite direction on the reverse strand.
data Range = Range {
        r_pos :: !Position, 
        r_length :: !Int 
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
shift :: Position -> Int -> Position
shift p a = case p_sense p of Forward -> p { p_start = p_start p + a }
                              Reverse -> p { p_start = p_start p - a }

-- | extends a range
-- The length of the range is simply increased.
extend :: Range -> Int -> Range
extend r a = r { r_length = r_length r + a }

