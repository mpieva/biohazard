{-# LANGUAGE RecordWildCards, NamedFieldPuns, BangPatterns, MultiParamTypeClasses, TypeFamilies, TemplateHaskell #-}
-- Estimates aDNA damage.  Crude first version.
--
-- - Read a BAM file, make compact representation of the reads.
-- - Compute likelihood of read under simple model of
--   damage, error/divergence, contamination.
--
-- For the less crude version:  follow the iobio subsampling strategy
-- using an index.

-- Trying to compute symbolically is too much, the high power terms get
-- out of hand quickly, and we get mixed powers of \lambda and \kappa.

-- For the fitting, we simplify radically: ignore sequencing error,
-- assume damage and simple, symmetric substitutions.

-- The fastest version so far uses the cheap implementation of automatic
-- differentiation in AD.hs together with the Hager-Zhang method from
-- nonlinear-optimization.  BFGS from hmatrix-gsl takes longer to
-- converge.
--
-- If I include parameters, whose true value is zero, the transformation
-- to the log-odds-ratio doesn't work, because then the maximum doesn't
-- exist anymore.  A different transformation ('sigmoid2'/'isigmoid2'
-- below) allows for an actual zero (but not one), while avoiding ugly
-- boundary conditions.  That works.
--
-- The current hack assumes all molecules have an overhang at both ends,
-- then each base gets deaminated with a position dependent probability.
-- If we try to model a fraction of undeaminated molecules in addition,
-- this fails.  To rescue the idea, I guess we must really decide if the
-- molecule has an overhang at all (probability 1/2) at each end, then
-- deaminate it.

-- TODO:
-- Before this can be packaged and used, the following needs to be done:
-- - Start with a crude estimate of parameters,
-- - Fix the model(s), so a contaminant fraction can be estimated.
-- - Implement both SSD and DSD.

-- final estimate [5.770143305791021e-25,0.4686003683523658,1.93782286620759e-2,0.6960685465309394,0.2780883283160801]


import Bio.Bam.Header
import Bio.Bam.Index
import Bio.Bam.Raw
import Bio.Bam.Rec
import Bio.Base
import Bio.Genocall.Adna
import Bio.Iteratee
import Control.Applicative
import Data.Bits
import Data.Foldable
import Data.Monoid
import Data.Traversable
import Data.Vec ( vec, Mat44, dot, getElem )
import Data.Vector.Unboxed.Deriving
import Numeric.Optimization.Algorithms.HagerZhang05
import System.Environment

import qualified Data.Vector                as V
import qualified Data.Vector.Fusion.Stream  as S
import qualified Data.Vector.Generic        as G
import qualified Data.Vector.Unboxed        as U

import AD
import Prelude hiding ( mapM_ )

import Debug.Trace

-- | Roughly @Maybe (Nucleotide, Nucleotide)@, encoded compactly
newtype NP = NP { unNP :: Word8 } deriving Eq
derivingUnbox "NP" [t| NP -> Word8 |] [| unNP |] [| NP |]

instance Show NP where
    show (NP w) | w > 15 = "NN"
    show (NP w) = [ "ACGT" !! fromIntegral (w `shiftR` 2)
                  , "ACGT" !! fromIntegral (w .&. 3) ]

isValid :: NP -> Bool
isValid (NP x) = x < 16

type Seq = U.Vector NP

sigmoid, sigmoid2, isigmoid, isigmoid2 :: (Num a, Fractional a, Floating a) => a -> a
sigmoid l = 1 / ( 1 + exp l )
isigmoid p = log $ ( 1 - p ) / p

sigmoid2 x = y*y where y = (exp x - 1) / (exp x + 1)
isigmoid2 y = log $ (1 + sqrt y) / (1 - sqrt y)

{-# INLINE lk_fun1 #-}
lk_fun1 :: (Num a, Fractional a, Floating a) => Seq -> [a] -> a
lk_fun1 bb parms = negate $ log $ lk bb
  where
    l_subst:l_sigma:l_delta:l_lam:l_kap:_ = parms

    -- Good initial guesses may be necessary, too.
    -- f_exo = sigmoid l_endo
    p_self = 1 - sigmoid2 l_subst
    p_subst = (1 - p_self) * 0.333

    ssd_sigma  = sigmoid l_sigma
    ssd_delta  = sigmoid l_delta
    ssd_lambda = sigmoid l_lam
    ssd_kappa  = sigmoid l_kap

    -- Technically, its transpose.  But it's symmetric anyway.
    subst_mat = vec4 (vec4 p_self p_subst p_subst p_subst)
                     (vec4 p_subst p_self p_subst p_subst)
                     (vec4 p_subst p_subst p_self p_subst)
                     (vec4 p_subst p_subst p_subst p_self)

    lk br = {-(1-f_exo) *-} S.foldl' (*) 1 (S.zipWith lk1 (G.stream br) $
                                          G.stream (ssDamage SSD{..} False (U.length br)))
            -- + f_exo * G.foldl' (\acc pr -> acc * lk1 pr W.identity) 1 br

    lk1 (NP pr) m | pr > 15   = 1
                  | otherwise = getElem b m `dot` getElem a subst_mat
      where
        a = fromIntegral $ pr `shiftR` 2  -- from
        b = fromIntegral $ pr .&. 3       -- to


lkfun :: V.Vector Seq -> U.Vector Double -> Double
lkfun brs parms = V.foldl' (\a b -> a + lk_fun1 b ps) 0 brs
  where
    !ps = U.toList parms

combofn :: V.Vector Seq -> U.Vector Double -> (Double, U.Vector Double)
combofn brs parms = (x,g)
  where
    !ps     = paramVector $ U.toList parms
    !n      = U.length parms
    (D x g) = V.foldl' (\a b -> a + lk_fun1 b ps) 0 brs


main :: IO ()
main = do
    [fp] <- getArgs
    brs <- subsampleBam fp >=> run $ \_ ->
           joinI $ filterStream (not . br_isUnmapped) $
           joinI $ mapStream pack_record $
           joinI $ filterStream (U.all isValid) $
           stream2vectorN 10000

    putStrLn $ "no. input sequences " ++ show (V.length brs)
    case crude_estimate brs of
        Left  v0 -> do putStrLn $ "crude estimate (SS): " ++ show (backxform v0)

                       let params = defaultParameters { verbose = VeryVerbose }
                       (xs, r, st) <- optimize params 0.00001 v0
                                               (VFunction $ lkfun brs)
                                               (VGradient $ snd . combofn brs)
                                               (Just . VCombined $ combofn brs)

                       print r
                       print st
                       putStrLn $ "crude estimate (SS): " ++ show (backxform v0)
                       putStrLn $ "final estimate (SS): " ++ show (backxform xs)

        Right v0 -> putStrLn $ "crude estimate (DS): " ++ show (backxform v0)

backxform vec = map sigmoid2 (G.toList $ G.take 1 vec)
             ++ map sigmoid (G.toList $ G.drop 1 vec)


-- We'll require the MD field to be present.  Then we cook each read
-- into a list of paired bases.  Deleted bases are dropped, inserted
-- bases replaced with an escape code.
--
-- XXX  This is annoying... almost, but not quite the same as the code
-- in the "Pileup" module.  This also relies on MD and doesn't offer the
-- alternative of accessing a reference genome.  (The latter may not be
-- worth the trouble.)  It also resembles the 'ECig' logic from
-- "Bio.Bam.Rmdup".

pack_record :: BamRaw -> Seq
pack_record br = if br_isReversed br then revcom u1 else u1
  where
    BamRec{..} = decodeBamEntry br

    revcom = U.reverse . U.map (NP . xor 15 . unNP)

    u1 = U.fromList $ go (unCigar b_cigar) (U.toList b_seq) (br_get_md br)
    go :: [(CigOp,Int)] -> [Nucleotides] -> [MdOp] -> [NP]

    go ((_,0):cs)   ns mds  = go cs ns mds
    go cs ns (MdNum  0:mds) = go cs ns mds
    go cs ns (MdDel []:mds) = go cs ns mds
    go  _ []              _ = []

    go ((Mat,nm):cs) (n:ns) (MdNum mm:mds) = mk_pair n n  : go ((Mat,nm-1):cs) ns (MdNum (mm-1):mds)
    go ((Mat,nm):cs) (n:ns) (MdRep n':mds) = mk_pair n n' : go ((Mat,nm-1):cs) ns               mds
    go ((Mat,nm):cs)    ns  (MdDel ds:mds) =                go ((Mat, nm ):cs) ns               mds

    go ((Ins,nm):cs) ns mds = replicate nm esc ++ go cs (drop nm ns) mds
    go ((SMa,nm):cs) ns mds = replicate nm esc ++ go cs (drop nm ns) mds
    go ((Del,nm):cs) ns (MdDel (_:ds):mds) = go ((Del,nm-1):cs) ns (MdDel ds:mds)
    go ((Del,nm):cs) ns (           _:mds) = go ((Del, nm ):cs) ns           mds

    go (_:cs) nd mds = go cs nd mds


esc :: NP
esc = NP 16

mk_pair :: Nucleotides -> Nucleotides -> NP
mk_pair (Ns a) = case a of 1 -> mk_pair' 0
                           2 -> mk_pair' 1
                           4 -> mk_pair' 2
                           8 -> mk_pair' 3
                           _ -> const esc
  where
    mk_pair' a (Ns b) = case b of 1 -> NP $ a .|. 0
                                  2 -> NP $ a .|. 4
                                  4 -> NP $ a .|. 8
                                  8 -> NP $ a .|. 12
                                  _ -> esc


infix 6 /%/
(/%/) :: Integral a => a -> a -> Double
0 /%/ 0 = 0
a /%/ b = fromIntegral a / fromIntegral (a+b)

-- Crude estimate.  Need two overhang lengths, two deamination rates,
-- undamaged fraction, SS/DS, substitution rate.
--
-- DS or SS: look whether CT or GA is greater at 3' terminal position  √
-- Left overhang length:  ratio of damage at second position to first  √
-- Right overang length:  ratio of CT at last to snd-to-last posn      √
--                      + ratio of GA at last to snd-to-last posn      √
-- SS rate: condition on damage on one end, compute rate at other      √
-- DS rate: condition on damage, compute rate in interior              √
-- substitution rate:  count all substitutions not due to damage       √
-- undamaged fraction:  ???


crude_estimate :: V.Vector Seq -> Either (U.Vector Double) (U.Vector Double)
crude_estimate seqs0
    | rate_ct > rate_ga = Left  $ U.fromList [ l_subst, l_sigma, l_delta, l_lam, l_kap ]        -- SS
    | otherwise         = Right $ U.fromList [ l_subst, l_sigma, l_delta, l_lam ]               -- DS
  where
    seqs = V.filter ((>= 10) . U.length) seqs0

    total_equals = V.sum (V.map (U.length . U.filter      isNotSubst) seqs)
    total_substs = V.sum (V.map (U.length . U.filter isOrdinarySubst) seqs) * 6 `div` 5
    l_subst = isigmoid2 $ total_substs /%/ total_equals

    isNotSubst (NP x) = x < 16 && x `shiftR` 2 == x .&. 3
    isOrdinarySubst (NP x) = x < 16 && x `shiftR` 2 /= x .&. 3 &&
                             x /= 7 {- CT -} && x /= 8 {- GA -}

    ct_at_alpha = V.length $ V.filter (\v -> v U.! 0 == NP 7 && dmg_omega v) seqs
    cc_at_alpha = V.length $ V.filter (\v -> v U.! 0 == NP 5 && dmg_omega v) seqs
    ct_at_beta  = V.length $ V.filter (\v -> v U.! 1 == NP 7 && dmg_omega v) seqs
    cc_at_beta  = V.length $ V.filter (\v -> v U.! 1 == NP 5 && dmg_omega v) seqs

    dmg_omega v = v U.! (l-1) == NP 7 || v U.! (l-1) == NP 8 ||
                  v U.! (l-2) == NP 7 || v U.! (l-2) == NP 8 ||
                  v U.! (l-3) == NP 7 || v U.! (l-3) == NP 8
        where l = U.length v

    l_lam = isigmoid lambda
    lambda = (ct_at_beta /%/ cc_at_beta) / (ct_at_alpha /%/ cc_at_alpha)

    ct_at_omega = V.length $ V.filter (\v -> v U.! (U.length v -1) == NP 7 && dmg_alpha v) seqs
    cc_at_omega = V.length $ V.filter (\v -> v U.! (U.length v -1) == NP 5 && dmg_alpha v) seqs
    ct_at_psi   = V.length $ V.filter (\v -> v U.! (U.length v -2) == NP 7 && dmg_alpha v) seqs
    cc_at_psi   = V.length $ V.filter (\v -> v U.! (U.length v -2) == NP 5 && dmg_alpha v) seqs

    ga_at_omega = V.length $ V.filter (\v -> v U.! (U.length v -1) == NP  8 && dmg_alpha v) seqs
    gg_at_omega = V.length $ V.filter (\v -> v U.! (U.length v -1) == NP 10 && dmg_alpha v) seqs

    dmg_alpha v = v U.! 0 == NP 7 || v U.! 1 == NP 7 || v U.! 2 == NP 7

    l_kap = isigmoid $ (ct_at_psi /%/ cc_at_psi) / (ct_at_omega /%/ cc_at_omega)

    rate_ct = ct_at_omega /%/ cc_at_omega
    rate_ga = ga_at_omega /%/ gg_at_omega

    total_inner_CCs = V.sum $ V.map (U.length . U.filter (== NP 5) . takeInner) seqs
    total_inner_CTs = V.sum $ V.map (U.length . U.filter (== NP 7) . takeInner) seqs
    takeInner v = U.slice 5 (U.length v - 10) v

    l_delta = isigmoid $ total_inner_CTs /%/ total_inner_CCs
    raw_rate = ct_at_alpha /%/ cc_at_alpha
    l_sigma = isigmoid $ {- XXX  2 * -} raw_rate / lambda

