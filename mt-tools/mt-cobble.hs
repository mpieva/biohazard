{-# LANGUAGE OverloadedStrings, BangPatterns #-}
{-# OPTIONS_GHC -Wall #-}

-- Cobble up a mitochondrion, or something similar.
--
-- The goal is to reconstruct a mitochondrion (or similar small, haploid
-- locus) from a set of sequencing reads and a reference sequence.  The
-- idea is to first select reads using some sort of filtering strategy,
-- simply for speed reasons.  They are then aligned to the reference
-- using banded Smith-Waterman algorithm, and a more likely reference is
-- called.  This is repeated till it converges.  A bad implementation of
-- the idea was called MIA.

import Align
import SimpleSeed

import Bio.Base
import Bio.Bam.Raw
import Bio.Iteratee
import Control.Applicative
import Data.Char
import System.Environment

import qualified Data.IntMap as IM
import qualified Data.ByteString.Lazy.Char8 as L


-- Read a FastA file, drop the names, yield the sequences.
readFasta :: L.ByteString -> [[Either Nucleotide Nucleotide]]
readFasta = go . L.lines
  where
    isHeader s = not (L.null s) && L.head s == '>'
    go ls = case break isHeader $ dropWhile isHeader ls of
                (body, rest) -> let ns = map toNuc . concat $ map L.unpack body
                                in ns : if null rest then [] else go rest
    toNuc x | isUpper x = Right $ toNucleotide x
            | otherwise = Left  $ toNucleotide (toUpper x)


main :: IO ()
main = do
    s:rs <- getArgs
    inputs@(reference:_) <- concatMap readFasta <$> mapM L.readFile rs
    let !ln = length reference
        !sm = create_seed_maps (map (map (either id id)) inputs)
        !rs = prep_reference reference

    print (ln, IM.size sm)
    -- print rs

    decodeAnyBamFile s >=> run $ \_ -> mapStreamM_ (round1 ln sm rs)


round1 ln sm rs br = case do_seed ln sm br of
    Nothing    -> return ()
    Just (a,b) | a >= 0 -> let qs = prep_query_fwd br
                           in test qs (RP a) (BW $ b-a-br_l_seq br)
               | otherwise -> let qs = revcompl_query $ prep_query_fwd br
                              in test qs (RP (-b)) (BW $ b-a-br_l_seq br)
  where
    test qs (RP x) (BW y) = do putStrLn $ "Rgn " ++ show x ++ ".." ++ show (x+br_l_seq br) ++ "x" ++ show y
                               print qs
                               let memo = viterbi_forward 50 rs qs (RP x) (BW y)
                               pmax (BW y) memo

