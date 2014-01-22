{-# LANGUAGE OverloadedStrings, BangPatterns #-}
{-# OPTIONS_GHC -Wall #-}

module Align where

import Bio.Base
import Bio.Bam
import Control.Monad
import Data.Bits
import Data.Int
import Data.Vector.Unboxed ((!))

import qualified Data.Vector.Unboxed         as U
import qualified Data.Vector.Unboxed.Mutable as UM

data Base = A | C | G | T | None
  deriving (Eq, Ord, Enum, Show)

-- | For a reference sequence, we store five(!) probabilities for each
-- base in phred format.  The fifth is the probability of a gap.

newtype RefSeq = RS (U.Vector Word8) deriving Show

refseq_len :: RefSeq -> Int
refseq_len (RS v) = U.length v `div` 5

prob_of :: Base -> Int -> RefSeq -> Word8
prob_of b i (RS v) = v ! ( 5*i + fromEnum b )

-- | Turns a sequence into probabilities.  @Right n@ is an ordinary
-- 'Nucleotide', @Left n@ is one we think might be absent (e.g. because
-- it was soft masked in the input).
prep_reference :: [Either Nucleotide Nucleotide] -> RefSeq
prep_reference = RS . U.concat .  map (either (to probG) (to probB))
  where
    to ps n = U.slice (5 * fromIntegral (unN n)) 5 ps

    probB = U.fromListN 80 $ concatMap (\l ->          l ++ [255]) raw_probs
    probG = U.fromListN 80 $ concatMap (\l -> map (+3) l ++  [3])  raw_probs

    raw_probs = [[ 252, 252, 252, 252 ]    -- 0
                ,[   0, 252, 252, 252 ]    -- A
                ,[ 252,   0, 252, 252 ]    -- C
                ,[   3,   3, 252, 252 ]    -- M
                ,[ 252, 225,   0, 252 ]    -- G
                ,[   3, 252,   3, 252 ]    -- R
                ,[ 252,   3,   3, 252 ]    -- S
                ,[   5,   5,   5, 252 ]    -- V
                ,[ 252, 252, 252,   0 ]    -- T
                ,[   3, 252, 252,   3 ]    -- W
                ,[ 252,   3, 252,   3 ]    -- Y
                ,[   5,   5, 252,   5 ]    -- H
                ,[ 252, 252,   3,   3 ]    -- K
                ,[   5, 252,   5,   5 ]    -- D
                ,[ 252,   5,   5,   5 ]    -- B
                ,[   6,   6,   6,   6 ]]   -- N

-- | Encoding of the query:  one word per position, the two lowest bits
-- encode the base, the rest is the quality score (shifted left by 2).
newtype QuerySeq = QS (U.Vector Word8) deriving Show

-- | Prepare query for subsequent alignment to the forward strand.
prep_query_fwd :: BamRaw -> QuerySeq
prep_query_fwd br = QS $ U.fromListN len
    [ q `shiftL` 2 .|. code ! b | i <- [0 .. len-1]
    , let b = fromIntegral $ unN $ br_seq_at  br i
    , let q = fromIntegral $ unQ $ br_qual_at br i ]
  where
    len  = br_l_seq br
    code = U.fromListN 16 [0,0,1,0,2,0,0,0,3,0,0,0,0,0,0,0]

prep_query_rev :: BamRaw -> QuerySeq
prep_query_rev = revcompl_query . prep_query_fwd
  where
  revcompl_query (QS v) = QS $ U.map (xor 3) $ U.reverse v


-- | Memoization matrix for dynamic programming.  We understand it as a
-- matrix B columns wide and L rows deep, where B is the bandwidth and L
-- the query length.  Successive rows are understood to be skewed to the
-- right.  (This means all operations need the bandwidth as an
-- argument.)

newtype MemoMat   = MemoMat (U.Vector Int32)
newtype Bandwidth = BW Int
newtype RefPosn   = RP Int

get_memo_max :: Bandwidth -> MemoMat -> Int
get_memo_max (BW bw) (MemoMat v) = fromIntegral $ U.minimum $ U.drop (U.length v - bw) v

-- | Smith-Waterman, banded, linear gap costs.
viterbi_forward :: Int32 -> RefSeq -> QuerySeq -> RefPosn -> Bandwidth -> MemoMat
viterbi_forward gp (RS rs) (QS qs) (RP p0) (BW bw) = MemoMat $ U.create (do
    v <- UM.unsafeNew $ bw * U.length qs

    let score qpos rpos | qpos < 0 || qpos >= U.length qs = error $ "Read from QS: " ++ show qpos
        score qpos rpos = let base = (qs ! qpos) .&. 3 :: Word8
                              qual = (qs ! qpos) `shiftR` 2 :: Word8
                              prob = let ix = 5*rpos + fromIntegral base in
                                     if ix < 0 then error ("Huh? " ++ show ix) else
                                     if ix < U.length rs then rs ! ix else
                                     if ix - U.length rs < U.length rs then rs ! (ix - U.length rs) :: Word8 else
                                     error ("Huh? " ++ show (ix,qpos,rpos,p0,base))
                          in fromIntegral $ min qual prob

    let gscore rpos = let prob = let ix = 5*rpos + 4 in
                                 if ix < 0 then error ("Huh? " ++ show ix) else
                                     if ix < U.length rs then rs ! ix else
                                     if ix - U.length rs < U.length rs then rs ! (ix - U.length rs) :: Word8 else
                                     error ("Huh? " ++ show (ix,rpos,p0))
                          in min gp $ fromIntegral prob

    let readV row col | row < 0 || col < 0 || col >= bw || row >= U.length qs = error $ "Read from memo: " ++ show (row,col)
                      | ix < 0 || ix >= UM.length v                           = error $ "Read from memo: " ++ show ix
                      | otherwise = UM.read v ix
            where ix = bw*row + col

    let match row col = (+ score (row-1) (p0+col)) `liftM` readV (row-1) (col+0)
        gapH  row col = (+ gscore (p0+col))        `liftM` readV (row+0) (col-1)
        gapV  row col = (+ gp)                     `liftM` readV (row-1) (col+1)

    -- first line
    forM_ [0] $ \row -> do
        -- first cell
        let col = 0 in
            UM.write v (bw*row + col) =<< return 0

        -- most cells
        forM_ [1 .. bw-1] $ \col ->
            UM.write v (bw*row + col) =<< liftM (min 0) (gapH row col)

    -- the other lines
    forM_ [1 .. U.length qs-1] $ \row -> do
        -- first cell
        let col = 0 in
            UM.write v (bw*row + col) =<< liftM2 min  (match row col) (gapV row col)

        -- most cells
        forM_ [1 .. bw-2] $ \col ->
            UM.write v (bw*row + col) =<< liftM3 min3 (match row col) (gapH row col) (gapV row col)

        -- last cell
        let col = bw-1 in
            UM.write v (bw*row + col) =<< liftM2 min  (match row col) (gapH row col)

    return v)
  where
    min3 a b c = min a $ min b c
