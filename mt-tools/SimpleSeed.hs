{-# LANGUAGE OverloadedStrings, BangPatterns #-}
{-# OPTIONS_GHC -Wall #-}
module SimpleSeed where

import Bio.Base
import Bio.Bam.Raw
import Bio.Iteratee

import Data.Array.Unboxed
import Data.Bits
import Data.List
import Data.Maybe
import System.IO

import qualified Data.IntMap as IM
import qualified Data.ByteString.Char8 as S

-- | Discontiguous template "12 of 16", stolen from MegaBLAST:
-- 1,110,110,110,110,111, with two bits per base gives 0xFCF3CF3F

template :: Int
template = 0xFCF3CF3F

create_seed_words :: [Nucleotide] -> [(Int, Int)]
create_seed_words = drop 32 . go 0x0 (-16) 0x0 0
  where
    go !accf !i !accr !ir s =
        (accf .&. template, i) : (accr .&. template, ir) : case s of
            [    ] -> []
            (n:ns) -> go (accf `shiftR` 2 .|. (codef ! n) `shiftL` 30) (i+1)
                         (accr `shiftL` 2 .|. (coder ! n)) (ir-1) ns

    -- These codes are chosen so that ambiguity codes result in zeroes.
    -- The seed word 0, which would otherwise be the low-complexity and
    -- useless poly-A, is later ignored.
    codef, coder :: UArray Nucleotide Int
    codef = listArray (N 0, N 15) [0,0,1,0,2,0,0,0,3,0,0,0,0,0,0,0]
    coder = listArray (N 0, N 15) [0,3,2,0,1,0,0,0,0,0,0,0,0,0,0,0]

-- Turns a list of seed words into a map.  Only the first entry is used,
-- duplicates are discarded silenty.

create_seed_map :: [Nucleotide] -> IM.IntMap Int
create_seed_map = cleanup . IM.fromListWith add . map (\(x,y) -> (x,(1::Int,y))) . create_seed_words . pad
  where pad ns = ns ++ take 15 ns
        add (x,i) (y,_) = (x+y,i)
        cleanup = IM.mapMaybe (\(n,j) -> if n < 8 then Just j else Nothing)

create_seed_maps :: [[Nucleotide]] -> IM.IntMap Int
create_seed_maps = IM.unionsWith const . map create_seed_map

-- | Actual seeding.  We take every hit and guesstimate an alignment
-- region from it (by adding the overhanging sequence parts and rounding
-- a bit up).  Regions are overlapped into larger ones, counting votes.
-- The region with the most votes is used as seed region.  (This will
-- occasionally result in a very long initial alignment.  We can afford
-- that.)
--
-- If we have PE data where only one read is seeded, we can either
-- discard the pair or align the second mate very expensively.  While
-- possible, that sounds rather expensive and should probably depend on
-- the quality of the first mates' alignment.  Generally, we may want to
-- check the quality of the initial alignment anyway.
--
-- Note that the overlapping logic right now doesn't take care of the
-- origin.  A simple fix would confuse both strand... and I haven't come
-- up with a complicated fix yet.  This should become a simple function,
-- once sufficiently debugged.
--
-- For proper overlapping:  We need to normalize each region to strictly
-- positive or strictly negative coordinates.  After sorting and
-- overlapping, we only need to check if the last region overlaps the
-- first---there can be only one such overlap per strand.  We should
-- probably discard overly long regions.

do_seed :: Int -> IM.IntMap Int -> BamRaw -> IO ()
do_seed ln sm br = do S.hPut stdout $ S.concat [ br_qname br, key, ":  ", S.pack (shows br_seq "\n") ]
                      mapM_ (\x -> hPutStrLn stdout $ "  " ++ show x) rgns
                      case rgns of
                           [         ] -> putStrLn "discard"
                           (a,b,_) : _ -> putStrLn $ "seed to " ++ shows a ".." ++ shows b " ("
                                                             ++ shows (b-a) "/" ++ shows (br_l_seq br) ")"

  where
    seeds = filter ((/= 0) . fst) $ filter ((/= template) . fst) $
            filter ((>= 0) . snd) $ create_seed_words br_seq

    br_seq = [ br_seq_at br i | i <- [0..br_l_seq br-1] ]

    key | br_isPaired br && br_isFirstMate  br = "/1"
        | br_isPaired br && br_isSecondMate br = "/2"
        | otherwise                            = "/m"

    more x = (x * 17) `div` 16

    rgns = sortBy (\(_,_,c) (_,_,z) -> compare z c) $ overlap $ sort $
               [ (j - more i, j + more (br_l_seq br-i), 1::Int)
               | (k,i) <- seeds, j <- maybeToList $ IM.lookup k sm ]

    overlap ( (x,y,n) : (i,j,m) : rs ) | i <= y = overlap ( (x,max y j,n+m) : rs )
    overlap ( (x,y,n) : rs ) = (x,y,n) : overlap rs
    overlap [] = []


