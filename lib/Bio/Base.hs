{-# LANGUAGE GeneralizedNewtypeDeriving, TypeFamilies, FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses, BangPatterns #-}
-- | Common data types used everywhere.  This module is a collection of
-- very basic "bioinformatics" data types that are simple, but don't
-- make sense to define over and over.

module Bio.Base(
    Nucleotide(..), Nucleotides(..),
    Qual(..), toQual, fromQual, fromQualRaised, probToQual,
    Prob(..), toProb, fromProb, qualToProb, pow,

    Word8,
    nucA, nucC, nucG, nucT,
    nucsA, nucsC, nucsG, nucsT, nucsN, gap,
    toNucleotide, toNucleotides,
    showNucleotide, showNucleotides,
    isGap,
    isBase,
    isProperBase,
    properBases,
    compl, compls,
    revcompl, revcompls,
    everything,

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

import Bio.Util             ( log1p )
import Data.Array.Unboxed
import Data.Bits
import Data.Char            ( isAlpha, isSpace, ord, toUpper )
import Data.Word            ( Word8 )
import Foreign.Storable     ( Storable(..) )
import Numeric              ( showFFloat )
import System.Directory     ( doesFileExist )
import System.FilePath      ( (</>), isAbsolute, splitSearchPath )
import System.Environment   ( getEnvironment )

import qualified Data.ByteString.Char8 as S
import qualified Data.ByteString.Lazy as L

import qualified Data.Vector.Generic         as VG
import qualified Data.Vector.Generic.Mutable as VM
import qualified Data.Vector.Unboxed         as VU

import Data.ByteString.Internal ( c2w, w2c )

-- | A nucleotide base.  We only represent A,C,G,T, and it is statically
-- enforced.

newtype Nucleotide = N { unN :: Word8 } deriving
    ( Eq, Ord, Ix, Storable
    , VG.Vector VU.Vector, VM.MVector VU.MVector, VU.Unbox )

instance Bounded Nucleotide where
    minBound = N 0
    maxBound = N 3

everything :: (Bounded a, Ix a) => [a]
everything = range (minBound, maxBound)

-- | A nucleotide base in an alignment.
-- Experience says we're dealing with Ns and gaps all the type, so
-- purity be damned, they are included as if they were real bases.
--
-- To allow @Nucleotides@s to be unpacked and incorparated into
-- containers, we choose to represent them the same way as the BAM file
-- format:  as a 4 bit wide field.  Gaps are encoded as 0 where they
-- make sense, N is 15.

newtype Nucleotides = Ns { unNs :: Word8 } deriving
    ( Eq, Ord, Ix, Storable
    , VG.Vector VU.Vector, VM.MVector VU.MVector, VU.Unbox )

instance Bounded Nucleotides where
    minBound = Ns  0
    maxBound = Ns 15

-- | Qualities are stored in deciban, also known as the Phred scale.  To
-- represent a value @p@, we store @-10 * log_10 p@.  Operations work
-- directly on the \"Phred\" value, as the name suggests.  The same goes
-- for the 'Ord' instance:  greater quality means higher \"Phred\"
-- score, meand lower error probability.
newtype Qual = Q { unQ :: Word8 } deriving
    ( Eq, Ord, Storable, Bounded, VG.Vector VU.Vector, VM.MVector VU.MVector, VU.Unbox )

instance Show Qual where
    showsPrec p (Q q) = (:) 'q' . showsPrec p q

toQual :: (Floating a, RealFrac a) => a -> Qual
toQual a = Q $ round (-10 * log a / log 10)

fromQual :: Qual -> Double
fromQual (Q q) = 10 ** (- fromIntegral q / 10)

fromQualRaised :: Double -> Qual -> Double
fromQualRaised k (Q q) = 10 ** (- k * fromIntegral q / 10)

-- | A positive 'Double' value stored in log domain.  We store the
-- natural logarithm (makes computation easier), but allow conversions
-- to the familiar \"Phred\" scale used for 'Qual' values.
newtype Prob = Pr { unPr :: Double } deriving
    ( Eq, Ord, Storable, VG.Vector VU.Vector, VM.MVector VU.MVector, VU.Unbox )

instance Show Prob where
    showsPrec _ (Pr p) = (:) 'q' . showFFloat (Just 1) q
      where q = - 10 * p / log 10

instance Num Prob where
    fromInteger a = Pr (log (fromInteger a))
    Pr x + Pr y = Pr $ if x >= y then x + log1p (  exp (y-x)) else y + log1p (exp (x-y))
    Pr x - Pr y = Pr $ if x >= y then x + log1p (- exp (y-x)) else error "no negative error probabilities"
    Pr a * Pr b = Pr $ a + b
    negate    _ = Pr $ error "no negative error probabilities"
    abs       x = x
    signum    _ = Pr 0

instance Fractional Prob where
    fromRational a = Pr (log (fromRational a))
    Pr a  /  Pr b = Pr (a - b)
    recip  (Pr a) = Pr (negate a)

infixr 8 `pow`
pow :: Prob -> Double -> Prob
pow (Pr a) e = Pr (a*e)


toProb :: Double -> Prob
toProb p = Pr (log p)

fromProb :: Prob -> Double
fromProb (Pr q) = exp q

qualToProb :: Qual -> Prob
qualToProb (Q q) = Pr (- log 10 * fromIntegral q / 10)

probToQual :: Prob -> Qual
probToQual (Pr p) = Q (round (- 10 * p / log 10))

nucA, nucC, nucG, nucT :: Nucleotide
nucA = N 0
nucC = N 1
nucG = N 2
nucT = N 3

gap, nucsA, nucsC, nucsG, nucsT, nucsN :: Nucleotides
gap   = Ns 0
nucsA = Ns 1
nucsC = Ns 2
nucsG = Ns 4
nucsT = Ns 8
nucsN = Ns 15


-- | Sequence identifiers are ASCII strings
-- Since we tend to store them for a while, we use strict byte strings.
-- If you get a lazy bytestring from somewhere, use 'shelve' to convert
-- it for storage.  Use @unpackSeqid@ and @packSeqid@ to avoid the
-- import of @Data.ByteString@.
type Seqid = S.ByteString

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


-- | Converts a character into a 'Nucleotides'.
-- The usual codes for A,C,G,T and U are understood, '-' and '.' become
-- gaps and everything else is an N.
toNucleotide :: Char -> Nucleotide
toNucleotide c = if inRange (bounds arr) (ord c) then N (arr ! ord c) else N 0
  where
    arr :: UArray Int Word8
    arr = listArray (0,127) (repeat 0) //
          ( [ (ord          x,  n) | (x, N n) <- pairs ] ++
            [ (ord (toUpper x), n) | (x, N n) <- pairs ] )

    pairs = [ ('a', nucA), ('c', nucC), ('g', nucG), ('t', nucT) ]


-- | Converts a character into a 'Nucleotides'.
-- The usual codes for A,C,G,T and U are understood, '-' and '.' become
-- gaps and everything else is an N.
toNucleotides :: Char -> Nucleotides
toNucleotides c = if inRange (bounds arr) (ord c) then Ns (arr ! ord c) else nucsN
  where
    arr :: UArray Int Word8
    arr = listArray (0,127) (repeat (unNs nucsN)) //
          ( [ (ord          x,  n) | (x, Ns n) <- pairs ] ++
            [ (ord (toUpper x), n) | (x, Ns n) <- pairs ] )

    Ns a `plus` Ns b = Ns (a .|. b)

    pairs = [ ('a', nucsA), ('c', nucsC), ('g', nucsG), ('t', nucsT),
              ('u', nucsT), ('-', gap),  ('.', gap),
              ('b', nucsC `plus` nucsG `plus` nucsT),
              ('d', nucsA `plus` nucsG `plus` nucsT),
              ('h', nucsA `plus` nucsC `plus` nucsT),
              ('v', nucsA `plus` nucsC `plus` nucsG),
              ('k', nucsG `plus` nucsT),
              ('m', nucsA `plus` nucsC),
              ('s', nucsC `plus` nucsG),
              ('w', nucsA `plus` nucsT),
              ('r', nucsA `plus` nucsG),
              ('y', nucsC `plus` nucsT) ]

-- | Tests if a 'Nucleotides' is a base.
-- Returns 'True' for everything but gaps.
isBase :: Nucleotides -> Bool
isBase (Ns n) = n /= 0

-- | Tests if a 'Nucleotides' is a proper base.
-- Returns 'True' for A,C,G,T only.
isProperBase :: Nucleotides -> Bool
isProperBase x = x == nucsA || x == nucsC || x == nucsG || x == nucsT

properBases :: [ Nucleotides ]
properBases = [ nucsA, nucsC, nucsG, nucsT ]

-- | Tests if a 'Nucleotides' is a gap.
-- Returns true only for the gap.
isGap :: Nucleotides -> Bool
isGap x = x == gap


{-# INLINE showNucleotide #-}
showNucleotide :: Nucleotide -> Char
showNucleotide (N x) = S.index str $ fromIntegral $ x .&. 3
  where str = S.pack "ACGT"

{-# INLINE showNucleotides #-}
showNucleotides :: Nucleotides -> Char
showNucleotides (Ns x) = S.index str $ fromIntegral $ x .&. 15
  where str = S.pack "-ACMGRSVTWYHKDBN"

instance Show Nucleotide where
    show     x = [ showNucleotide x ]
    showList l = (map showNucleotide l ++)

instance Read Nucleotide where
    readsPrec _ ('a':cs) = [(nucA, cs)]
    readsPrec _ ('A':cs) = [(nucA, cs)]
    readsPrec _ ('c':cs) = [(nucC, cs)]
    readsPrec _ ('C':cs) = [(nucC, cs)]
    readsPrec _ ('g':cs) = [(nucG, cs)]
    readsPrec _ ('G':cs) = [(nucG, cs)]
    readsPrec _ ('t':cs) = [(nucT, cs)]
    readsPrec _ ('T':cs) = [(nucT, cs)]
    readsPrec _ ('u':cs) = [(nucT, cs)]
    readsPrec _ ('U':cs) = [(nucT, cs)]
    readsPrec _     _    = [          ]

    readList ('-':cs) = readList cs
    readList (c:cs) | isSpace c = readList cs
                    | otherwise = case readsPrec 0 (c:cs) of
                            [] -> [ ([],c:cs) ]
                            xs -> [ (n:ns,r2) | (n,r1) <- xs, (ns,r2) <- readList r1 ]
    readList [] = [([],[])]

instance Show Nucleotides where
    show     x = [ showNucleotides x ]
    showList l = (map showNucleotides l ++)

instance Read Nucleotides where
    readsPrec _ (c:cs) = [(toNucleotides c, cs)]
    readsPrec _ [    ] = []
    readList s = let (hd,tl) = span (\c -> isAlpha c || isSpace c || '-' == c) s
                 in [(map toNucleotides $ filter (not . isSpace) hd, tl)]

-- | Reverse-complements a stretch of Nucleotides
{-# INLINE revcompl #-}
revcompl :: [Nucleotide] -> [Nucleotide]
revcompl = reverse . map compl

-- | Reverse-complements a stretch of Nucleotidess
{-# INLINE revcompls #-}
revcompls :: [Nucleotides] -> [Nucleotides]
revcompls = reverse . map compls

-- | Complements a Nucleotides.
{-# INLINE compl #-}
compl :: Nucleotide -> Nucleotide
compl (N n) = N $ n `xor` 3

-- | Complements a Nucleotides.
{-# INLINE compls #-}
compls :: Nucleotides -> Nucleotides
compls (Ns x) = Ns $ arr ! (x .&. 15)
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

