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
import Data.Sequence ( (<|), (><), ViewL((:<)) )

import qualified Data.ByteString             as S
import qualified Data.Foldable               as F
import qualified Data.Sequence               as Z
import qualified Data.Vector.Generic         as V
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
prob_of b i (RS v) = indexV "prob_of" v ( 5*i + fromEnum b )

-- | Turns a sequence into probabilities.  @Right n@ is an ordinary
-- 'Nucleotide', @Left n@ is one we think might be absent (e.g. because
-- it was soft masked in the input).
prep_reference :: [Either Nucleotides Nucleotides] -> RefSeq
prep_reference = RS . U.concat .  map (either (to probG) (to probB))
  where
    to ps n = U.slice (5 * fromIntegral (unNs n)) 5 ps

    -- XXX we should probably add some noise here, so the placement of
    -- gaps isn't completely random, but merely unpredictable
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
prep_query_fwd br = QS $ U.fromListN len $ zipWith pair (V.toList b_seq) (S.unpack b_qual)
  where
    BamRec{..} = unpackBam br
    pair b q = q `shiftL` 2 .|. indexV "prep_query_fwd" code (fromIntegral $ unNs b)
    code = U.fromListN 16 [0,0,1,0,2,0,0,0,3,0,0,0,0,0,0,0]
    len  = V.length b_seq

prep_query_rev :: BamRaw -> QuerySeq
prep_query_rev = revcompl_query . prep_query_fwd
  where
  revcompl_query (QS v) = QS $ U.map (xor 3) $ U.reverse v

qseqToBamSeq :: QuerySeq -> Vector_Nucs_half Nucleotides
qseqToBamSeq = V.fromList . U.toList . U.map (\x -> Ns $ 1 `shiftL` fromIntegral (x .&. 3)) . unQS

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
        { viterbi_forward :: MemoMat         -- DP matrix from running Viterbi
        , viterbi_score :: Float             -- alignment score (log scale, vs. radom alignment)
        , viterbi_position :: Int            -- position (start of the most probable alignment)
        , viterbi_backtrace :: Cigar }       -- backtrace (most probable alignment)
  deriving Show

data Traced = Tr { tr_op :: CigOp, tr_score :: Float }

instance Eq Traced where Tr _ a == Tr _ b = a == b
instance Ord Traced where Tr _ a `compare` Tr _ b = a `compare` b

-- | All sorts of alignment shit collected in one place, mostly so I can
-- reuse the scoring functions.
align :: Float -> RefSeq -> QuerySeq -> RefPosn -> Bandwidth -> AlignResult
align gp (RS rs) (QS qs) (RP p0) (BW bw_) = runST (do
    let bw = abs bw_
    v <- UM.unsafeNew $ bw * U.length qs + bw

    let readV row col | row < 0 || col < 0 || col >= bw || row > U.length qs = error $ "Read from memo: " ++ show (row,col)
                      | ix < 0 || ix >= UM.length v                          = error $ "Read from memo: " ++ show ix
                      | otherwise = UM.read v ix
            where ix = bw*row + col

    let score qpos    _ | qpos < 0 || qpos >= U.length qs = error $ "Read from QS: " ++ show qpos
        score qpos rpos = let base = (indexV "align/score/base" qs qpos) .&. 3 :: Word8
                              qual = (indexV "align/score/qual" qs qpos) `shiftR` 2 :: Word8
                              prob = let ix = 5*rpos + fromIntegral base in
                                     if ix < 0 then error ("Huh? " ++ show ix) else
                                     if ix < U.length rs then indexV "align/score/prob/A" rs ix else
                                     if ix - U.length rs < U.length rs then indexV "align/score/prob/B" rs (ix - U.length rs) :: Word8 else
                                     error ("Huh? " ++ show (ix,qpos,rpos,p0,base))

                              -- Improbability of a mismatch, it's the
                              -- probability of the reference not being
                              -- correct or the query not being correct,
                              -- whichever is higher.
                              mismatch = fromIntegral (min qual prob)

                              -- Improbability of a random match.  It's
                              -- 6 if we have a good base, corresponding
                              -- to randomness.  If we have a bad base,
                              -- it's lower, because we aren't doing
                              -- better than random.
                              randmatch = fromIntegral (min qual 6)

                              -- Score is our mismatch probability vs.
                              -- random sequences.  Note that this ends
                              -- up being 0 for low quality bases, -6
                              -- for high quality matches, and 30+ for
                              -- high quality mismatches.
                          in mismatch - randmatch

    let gscore rpos = let prob = let ix = 5*rpos + 4 in
                                 if ix < 0 then error ("Huh? " ++ show ix) else
                                     if ix < U.length rs then indexV "align/gscore/prob/A" rs ix else
                                     if ix - U.length rs < U.length rs then indexV "align/gscore/prob/B" rs (ix - U.length rs) :: Word8 else
                                     error ("Huh? " ++ show (ix,rpos,p0))
                          in min gp $ fromIntegral prob

    let match row col = Tr Mat . (+ score (row-1) (p0+row+col-1)) <$> readV (row-1) (col+0)
    let gapH  row col = Tr Del . (+ gscore (p0+row+col-1))        <$> readV (row+0) (col-1)
    let gapV  row col = Tr Ins . (+ gp)                           <$> readV (row-1) (col+1)

    let cell row col = do x <- if row == 0       then return (Tr Nop 0) else          match row col
                          y <- if             col == 0    then return x else min x <$> gapH row col
                          z <- if row == 0 || col == bw-1 then return y else min y <$> gapV row col
                          return z

    -- Fill the DP matrix.  XXX:  there's got to be way to express this
    -- using 'Vector's bulk operations.  Would that be more readable?
    -- Faster?
    forM_ [0 .. U.length qs] $ \row ->
        forM_ [0 .. bw-1] $ \col ->
            UM.write v (bw*row + col) . tr_score =<< cell row col

    let pack_cigar = Cigar . map (\x -> (head x, length x)) . group
    let traceback acc row col = do op <- tr_op <$> cell row col
                                   case op of Mat -> traceback (Mat:acc) (row-1) (col+0)
                                              Ins -> traceback (Ins:acc) (row-1) (col+1)
                                              Del -> traceback (Del:acc) (row+0) (col-1)
                                              Nop | row == 0 -> return (p0+col, pack_cigar acc)

    viterbi_forward <- MemoMat <$> U.unsafeFreeze v
    (viterbi_score, mincol) <- minimum . flip zip [0..] <$> mapM (readV (U.length qs)) [0..bw-1]
    (viterbi_position, viterbi_backtrace) <- traceback [] (U.length qs) mincol
    return $ AlignResult{..})

-- For each position, a vector of pseudocounts in the same order as in
-- 'RefSeq', followed by the same for based inserted after the current
-- one.
newtype NewRefSeq = NRS (Z.Seq NewColumn)

-- Inserts come (conceptually) before the base whose coordinate they
-- bear.  So every column has inserts first, then the single aligned
-- base.
data NewColumn = NC { nc_inserts :: !(U.Vector Float)
                    , nc_base    :: !(U.Vector Float) }

new_ref_seq :: RefSeq -> NewRefSeq
new_ref_seq rs = NRS $ Z.replicate (refseq_len rs) (NC (U.replicate 0 0) (U.replicate 5 0))

mkNC :: U.Vector Float -> U.Vector Float -> NewColumn
mkNC !i !b | U.length b /= 5 = error "mkNC"
           | otherwise = NC i b

-- Add an alignment to the new reference.  We compute the quality of the
-- alignment (probability that it belongs vs. probability that it's
-- random), that's how many votes we're going to cast.  (A perfect
-- alignment gives a whole vote, a random one gives none.  Call this
-- with an alignment that's worse than random at your own peril.)
-- If we're voting for a base, we vote for the called one according to
-- its quality and for all others with the error probability.
-- A deletion is a vote against all bases, an insert is a vote for how
-- ever many bases.  The first five values sum up to the total votes so
-- far, and they all count as votes against any further extension to an
-- insert.  We start with five pseudo-votes to get the numerics under
-- control (or to have a uniform Dirichlet-prior, if you prefer).
--
-- Note that this logic was arrived at by "thinking hard".  A clean way
-- to do it is to maximize the alignment score expected in the next
-- round, assuming the alignments do not change.  It might work out to
-- the same thing... who knows?

add_to_refseq :: NewRefSeq -> QuerySeq -> AlignResult -> NewRefSeq
add_to_refseq (NRS nrs0) (QS qs0) AlignResult{..} =
    NRS $ rotateZ (Z.length nrs0 - viterbi_position)
        $ mat here back qs0 (unCigar viterbi_backtrace)
  where
    here :< back = Z.viewl $ rotateZ viterbi_position nrs0
    rotateZ n = uncurry (flip (><)) . Z.splitAt n

    !odds = 10 ** (-viterbi_score / 10)  -- often huge,
    !votes = 1 - recip (1+odds)          -- often exactly 1

    -- Grrr, this isn't going to work.  We'll split it:
    -- One function deals with inserts.  As long as we get inserted
    -- bases, we vote for them.  Then we vote against the remainder and
    -- pass the buck.
    -- The other deals with a base.  We vote for it if we matched it,
    -- against it if we deleted it.  Then we recurse.
    ins !nc@(NC is b) !nrs !nins !qs cigs = case cigs of
        [          ] -> nc <| nrs
        (( _ ,0):cs) -> ins nc nrs nins qs cs

        ((Ins,n):cs) -> let is' = vote_for_at votes (U.sum b) nins (U.head qs) is
                        in ins (mkNC is' b) nrs (nins+1) (U.tail qs) ((Ins,n-1):cs)

        _            -> let is' = vote_against_from votes nins is
                        in mat (mkNC is' b) nrs qs cigs

    mat !nc@(NC is b) !nrs !qs cigs = case cigs of
        [          ] -> nc <| nrs
        (( _ ,0):cs) -> mat nc nrs qs cs

        ((Del,n):cs) -> let nc2 :< rest = Z.viewl nrs
                            b' = vote_against votes b
                        in mkNC is b' <!| mat nc2 rest qs ((Del,n-1):cs)

        ((Mat,n):cs) -> let nc2 :< rest = Z.viewl nrs
                            b' = vote_for votes (U.head qs) b
                        in mkNC is b' <!| mat nc2 rest (U.tail qs) ((Mat,n-1):cs)

        _            -> ins nc nrs (0::Int) qs cigs

    (<!|) !a !as = a <| as


vote_against_from :: Float -> Int -> U.Vector Float -> U.Vector Float
-- vote_against_from votes ix ps | trace ("vote_against_from " ++ show (ix, U.length ps)) False = undefined
vote_against_from votes ix ps = U.accum (+) ps [(i,votes) | i <- [ix+4, ix+9 .. U.length ps-1]]

vote_against :: Float -> U.Vector Float -> U.Vector Float
-- vote_against votes ps | trace ("vote_against " ++ show (U.length ps)) False = undefined
vote_against votes ps = U.accum (+) ps [(4,votes)]

vote_for :: Float -> Word8 -> U.Vector Float -> U.Vector Float
vote_for votes = vote_for_at votes 0 0

vote_for_at :: Float -> Float -> Int -> Word8 -> U.Vector Float -> U.Vector Float
vote_for_at votes v0 idx bq ps =
    U.accum (+) ps' $ (base+5*idx,pt) : [(5*idx+i,pe)|i<-[0,1,2,3]]
  where
    base = fromIntegral $ bq .&. 3
    qual = bq `shiftR` 2
    perr = 10 ** (fromIntegral qual * (-0.1))
    pe = votes * perr / 3
    pt = votes * (1 - perr) - pe

    ps' | U.length ps >= 5*idx+5 = ps
        | otherwise              = U.concat (ps : replicate (idx+1 - U.length ps `div` 5) (U.fromList [0,0,0,0,v0]))


-- Back to compact representation.  Every group of five votes gets
-- converted to five probabilities, and those to quality scores.  Then
-- we concatenate.
finalize_ref_seq :: NewRefSeq -> (RefSeq, XTab)
finalize_ref_seq (NRS z) =
    ( RS $ U.concat $ F.foldr unpack [] z
    , Z.fromList $ scanl (+) 0 $ F.foldr tolen [] z)
  where
    unpack (NC ins bas) k = map5 call ins ++ call bas : k
    map5 f v = [ f (U.slice i 5 v) | i <- [0, 5 .. U.length v - 5] ]
    call v = U.map (\x -> round $ (-10) / log 10 * log ((x+1) / total)) v where total = U.sum v + 5

    tolen (NC ins _bas) k = U.length ins `div` 5 + 1 : k

-- Table for coordinate translation
type XTab = Z.Seq Int


{-# INLINE indexV #-}
indexV :: String -> U.Vector Word8 -> Int -> Word8
-- indexV m v i | i  <          0 = error $ m ++ ": index too large"
             -- | i >= U.length v = error $ m ++ ": negative index"
             -- | otherwise       = v U.! i
indexV _ = (U.!)
