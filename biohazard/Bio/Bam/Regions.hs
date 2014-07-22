module Bio.Bam.Regions where

import Bio.Bam.Header ( Refseq(..) )
import Data.List ( foldl' )
import qualified Data.IntMap as IM

data Region = Region { rgn_refseq :: !Refseq, rgn_start :: !Int, rgn_end :: !Int }
  deriving (Eq, Ord, Show)

-- | A subset of a genome.  The idea is to map the reference sequence
-- (represented by its number) to a 'Subseqeunce'.
newtype Regions = Regions (IM.IntMap Subsequence) deriving Show

-- | A mostly contiguous subset of a sequence, stored as a set of
-- non-overlapping intervals in an 'IntMap' from start position to end
-- position (half-open intervals, naturally).
newtype Subsequence = Subsequence (IM.IntMap Int) deriving Show

toList :: Regions -> [(Refseq, Subsequence)]
toList (Regions m) = [ (Refseq $ fromIntegral k, v) | (k,v) <- IM.toList m ]

fromList :: [Region] -> Regions
fromList rs = Regions $ foldl' add_int IM.empty
    [ (fromIntegral r, beg, end) | Region (Refseq r) beg end <- rs ]
  where
    add_int m (r,b,e) = IM.alter (maybe single (Just . addInt b e)) r m
      where
        single = Just . Subsequence $ IM.singleton b e


addInt :: Int -> Int -> Subsequence -> Subsequence
addInt b e (Subsequence m0) = Subsequence $ merge_into b e m0
  where
    merge_into x y m = case IM.lookupLT y m of
        Just (u,v) | x < u && y <= v -> merge_into x v $ IM.delete u m    -- extend to the left
                   | x < u           -> merge_into x y $ IM.delete u m    -- subsume
                   | y <= v          -> m                                 -- subsumed
                   | x <= v          -> merge_into u y $ IM.delete u m    -- extend to the right
        _                            -> IM.insert  x y m                  -- no overlap

overlaps :: Int -> Int -> Subsequence -> Bool
overlaps b e (Subsequence m) = case IM.lookupLT e m of
        Just (_,v) -> b < v
        Nothing    -> False

