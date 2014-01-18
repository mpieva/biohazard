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
import System.Environment

import qualified Data.IntMap as IM
import qualified Data.ByteString.Lazy.Char8 as L


-- Read a FastA file, drop the names, yield the sequences.
readFasta :: L.ByteString -> [[Nucleotide]]
readFasta = go . L.lines
  where
    isHeader s = not (L.null s) && L.head s == '>'
    go ls = case break isHeader $ dropWhile isHeader ls of
                (body, rest) -> let ns = map toNucleotide . concat $ map L.unpack body
                                in ns : if null rest then [] else go rest


main :: IO ()
main = do
    s:rs <- getArgs
    inputs@(reference:_) <- concatMap readFasta <$> mapM L.readFile rs
    let !rl = length reference
        !sm = create_seed_maps inputs

    print (length i1, IM.size sm)

    decodeAnyBamFile s >=> run $ \_ -> mapStreamM_ (do_seed rl sm)


