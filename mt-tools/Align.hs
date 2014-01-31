{-# LANGUAGE OverloadedStrings, BangPatterns, RecordWildCards #-}
{-# OPTIONS_GHC -Wall #-}

module Align where

import Bio.Base
import Bio.Bam
import Control.Applicative
import Control.Monad
import Control.Monad.ST (runST)
import Data.Bits
import Data.List (group)
import Data.Vector.Unboxed ((!))

import qualified Data.ByteString             as S
import qualified Data.Vector.Unboxed         as U
import qualified Data.Vector.Unboxed.Mutable as UM

import Debug.Trace

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

    raw_probs = [[ 25, 25, 25, 25 ]    -- 0
                ,[  0, 25, 25, 25 ]    -- A
                ,[ 25,  0, 25, 25 ]    -- C
                ,[  3,  3, 25, 25 ]    -- M
                ,[ 25, 25,  0, 25 ]    -- G
                ,[  3, 25,  3, 25 ]    -- R
                ,[ 25,  3,  3, 25 ]    -- S
                ,[  5,  5,  5, 25 ]    -- V
                ,[ 25, 25, 25,  0 ]    -- T
                ,[  3, 25, 25,  3 ]    -- W
                ,[ 25,  3, 25,  3 ]    -- Y
                ,[  5,  5, 25,  5 ]    -- H
                ,[ 25, 25,  3,  3 ]    -- K
                ,[  5, 25,  5,  5 ]    -- D
                ,[ 25,  5,  5,  5 ]    -- B
                ,[  6,  6,  6,  6 ]]   -- N

-- | Encoding of the query:  one word per position, the two lowest bits
-- encode the base, the rest is the quality score (shifted left by 2).
newtype QuerySeq = QS { unQS :: U.Vector Word8 } deriving Show

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

qseqToBamSeq :: QuerySeq -> Sequence
qseqToBamSeq = U.map (\x -> N $ 1 `shiftL` fromIntegral (x .&. 3)) . unQS

qseqToBamQual :: QuerySeq -> S.ByteString
qseqToBamQual = S.pack . U.toList . U.map (`shiftR` 2) . unQS

-- | Memoization matrix for dynamic programming.  We understand it as a
-- matrix B columns wide and L rows deep, where B is the bandwidth and L
-- the query length.  Successive rows are understood to be skewed to the
-- right.  (This means all operations need the bandwidth as an
-- argument.)

newtype MemoMat   = MemoMat (U.Vector Float) deriving Show
newtype Bandwidth = BW Int deriving Show
newtype RefPosn   = RP Int deriving Show

data AlignResult = AlignResult
        { viterbi_forward :: MemoMat                    -- DP matrix from running Viterbi
        , viterbi_score :: Float
        , viterbi_backtrace :: (Int, Cigar) }       -- backtrace (most probable alignment)
  deriving Show

data Traced = Tr { tr_op :: CigOp, tr_score :: Float }

instance Eq Traced where Tr _ a == Tr _ b = a == b
instance Ord Traced where Tr _ a `compare` Tr _ b = a `compare` b

-- | All sorts of alignment shit collected in one place, mostly so I can
-- reuse the scoring functions.
align :: Float -> RefSeq -> QuerySeq -> RefPosn -> Bandwidth -> AlignResult
align gp (RS rs) (QS qs) (RP p0) (BW bw) = runST (do
    v <- UM.unsafeNew $ bw * U.length qs + bw

    let readV row col | row < 0 || col < 0 || col >= bw || row > U.length qs = error $ "Read from memo: " ++ show (row,col)
                      | ix < 0 || ix >= UM.length v                          = error $ "Read from memo: " ++ show ix
                      | otherwise = UM.read v ix
            where ix = bw*row + col

    let score qpos rpos | qpos < 0 || qpos >= U.length qs = error $ "Read from QS: " ++ show qpos
        score qpos rpos = let base = (qs ! qpos) .&. 3 :: Word8
                              qual = (qs ! qpos) `shiftR` 2 :: Word8
                              prob = let ix = 5*rpos + fromIntegral base in
                                     if ix < 0 then error ("Huh? " ++ show ix) else
                                     if ix < U.length rs then rs ! ix else
                                     if ix - U.length rs < U.length rs then rs ! (ix - U.length rs) :: Word8 else
                                     error ("Huh? " ++ show (ix,qpos,rpos,p0,base))
                          in fromIntegral (min qual prob) - 6  -- correcting for random matches

    let gscore rpos = let prob = let ix = 5*rpos + 4 in
                                 if ix < 0 then error ("Huh? " ++ show ix) else
                                     if ix < U.length rs then rs ! ix else
                                     if ix - U.length rs < U.length rs then rs ! (ix - U.length rs) :: Word8 else
                                     error ("Huh? " ++ show (ix,rpos,p0))
                          in min gp $ fromIntegral prob

    let match row col = Tr Mat . (+ score (row-1) (p0+row+col-1)) <$> readV (row-1) (col+0)
    let gapH  row col = Tr Del . (+ gscore (p0+row+col-1))        <$> readV (row+0) (col-1)
    let gapV  row col = Tr Ins . (+ gp)                           <$> readV (row-1) (col+1)

    let cell row col = do x <- if row == 0       then return (Tr Nop 0) else          match row col
                          y <- if             col == 0    then return x else min x <$> gapH row col
                          z <- if row == 0 || col == bw-1 then return y else min y <$> gapV row col
                          return z

    -- fill the DP matrix
    forM_ [0 .. U.length qs] $ \row ->
        forM_ [0 .. bw-1] $ \col ->
            UM.write v (bw*row + col) . tr_score =<< cell row col

    let pack_cigar = Cigar . map (\x -> (head x, length x)) . group
    let trace acc row col = do op <- tr_op <$> cell row col
                               case op of Mat -> trace (Mat:acc) (row-1) (col+0)
                                          Ins -> trace (Ins:acc) (row-1) (col+1)
                                          Del -> trace (Del:acc) (row+0) (col-1)
                                          Nop | row == 0 -> return (p0+col, pack_cigar acc)

    viterbi_forward <- MemoMat <$> U.unsafeFreeze v
    (viterbi_score, mincol) <- minimum . flip zip [0..] <$> mapM (readV (U.length qs)) [0..bw-1]
    viterbi_backtrace <- trace [] (U.length qs) mincol
    return $ AlignResult{..})
