{-# LANGUAGE OverloadedStrings, BangPatterns #-}
{-# OPTIONS_GHC -Wall #-}
module SimpleSeed where

import Bio.Base
import Bio.Bam.Raw

import Data.Array.Unboxed
import Data.Bits
import Data.List
import Data.Maybe

import qualified Data.IntMap as IM

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

data I2 = I2 !Int !Int

newtype SeedMap = SM { unSM :: IM.IntMap Int }
  deriving Show

create_seed_map ::  [Nucleotide] -> SeedMap
create_seed_map = SM . cleanup . foldl' (\m (k,v) -> IM.insertWith' add k v m) IM.empty .
                  map (\(x,y) -> (x,(I2 1 y))) . create_seed_words . pad
  where pad ns = ns ++ take 15 ns
        add (I2 x i) (I2 y _) = I2 (x+y) i
        cleanup = IM.mapMaybe $ \(I2 n j) -> if n < 8 then Just j else Nothing

create_seed_maps :: [[Nucleotide]] -> SeedMap
create_seed_maps = SM . IM.unionsWith const . map (unSM . create_seed_map)

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
-- For proper overlapping, we need to normalize each region to strictly
-- positive or strictly negative coordinates.  After sorting and
-- overlapping, we only need to check if the last region overlaps the
-- first---there can be only one such overlap per strand.  We should
-- probably discard overly long regions.

do_seed :: Int -> SeedMap -> BamRaw -> Maybe (Int,Int)
do_seed ln (SM sm) br = -- do S.hPut stdout $ S.concat [ br_qname br, key, ":  ", S.pack (shows br_seq "\n") ]
                   --    mapM_ (\x -> hPutStrLn stdout $ "  " ++ show x) rgns
                   case rgns of
                           [         ] -> Nothing -- putStrLn "discard"
                           {- (a,b,_) : _ | a > 20000 || a < (-20000) -> error $ concat [
                                    "Weird region: ",
                                    shows (a,b,ln) "; ",
                                    "Primitive regions: ",
                                    shows (rgns_fwd ++ rgns_rev) "; ",
                                    "Resulting regions: ",
                                    show rgns ] -}
                           (a,b,_) : _ -> Just (a,b) -- putStrLn $ "seed to " ++ shows a ".." ++ shows b " ("
                                                         --     ++ shows (b-a) "/" ++ shows (br_l_seq br) ")"

  where
    seeds = filter ((/= 0) . fst) $ filter ((/= template) . fst) $
            filter ((>= 0) . snd) $ create_seed_words br_seq

    br_seq = [ br_seq_at br i | i <- [0..br_l_seq br-1] ]

    more x = (x * 9) `div` 8 + 16

    rgns = sortBy (\(_,_,c) (_,_,z) -> compare z c) $ filter reasonably_short $
                (wrap_with        id $ overlap $ sort $ map norm_right rgns_fwd) ++
                (wrap_with norm_left $ overlap $ sort $ map norm_left  rgns_rev)

    (rgns_fwd, rgns_rev) = let put (f,r) (i,j) | j >= 0    = (rgn:f, r)
                                               | otherwise = (f, rgn:r)
                                where rgn = (j - more i, j + more (br_l_seq br-i), 1::Int)
                           in foldl put ([],[]) [ (i,j) | (k,i) <- seeds, j <- maybeToList $ IM.lookup k sm ]

    norm_right (a,b,n) = if a  < 0 then (a+ln, b+ln, n) else (a,b,n)
    norm_left  (a,b,n) = if b >= 0 then (a-ln, b-ln, n) else (a,b,n)

    wrap_with _    [           ] = []
    wrap_with _    [     r     ] = [r]
    wrap_with f rs@((x,y,n):rs')
        | i <= y+ln && x+ln <= j = f (min (x+ln) i, max (y+ln) j, n+m) : init rs'
        | otherwise              = rs
      where
        (i,j,m) = last rs

    overlap ( (x,y,n) : (i,j,m) : rs ) | i <= y = overlap ( (x,max y j,n+m) : rs )
    overlap ( (x,y,n) : rs ) = (x,y,n) : overlap rs
    overlap [] = []

    -- First cut:  reasonable is less than the whole MT.  Tuning can
    -- come later.
    reasonably_short (x,y,_) = y-x < ln

