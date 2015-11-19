{-# LANGUAGE RecordWildCards, NamedFieldPuns, BangPatterns, TypeFamilies #-}
-- Estimates aDNA damage.  Crude first version.
--
-- - Read or subsample a BAM file, make compact representation of the reads.
-- - Compute likelihood of each read under simple model of
--   damage, error/divergence, contamination.
--
-- For the fitting, we simplify radically: ignore sequencing error,
-- assume damage and simple, symmetric substitutions which subsume error
-- and divergence.
--
-- Trying to compute symbolically is too much, the high power terms get
-- out of hand quickly, and we get mixed powers of \lambda and \kappa.
-- The fastest version so far uses the cheap implementation of automatic
-- differentiation in AD.hs together with the Hager-Zhang method from
-- package nonlinear-optimization.  BFGS from hmatrix-gsl takes longer
-- to converge.  Didn't try an actual Newton iteration (yet?), AD from
-- package ad appears slower.
--
-- If I include parameters, whose true value is zero, the transformation
-- to the log-odds-ratio doesn't work, because then the maximum doesn't
-- exist anymore.  For many parameters, zero makes sense, but one
-- doesn't.  A different transformation ('sigmoid2'/'isigmoid2'
-- below) allows for an actual zero (but not an actual one), while
-- avoiding ugly boundary conditions.  That appears to work well.
--
-- The current hack assumes all molecules have an overhang at both ends,
-- then each base gets deaminated with a position dependent probability
-- following a geometric distribution.  If we try to model a fraction of
-- undeaminated molecules (a contaminant) in addition, this fails.  To
-- rescue the idea, I guess we must really decide if the molecule has an
-- overhang at all (probability 1/2) at each end, then deaminate it.
--
-- TODO
--   - needs better packaging, better output
--   - needs support for multiple input files(?)
--   - needs read group awareness(?)
--   - needs to deal with long (unmerged) reads (by ignoring them?)

import Bio.Bam.Header
import Bio.Bam.Index
import Bio.Bam.Raw
import Bio.Bam.Rec
import Bio.Base
import Bio.Genocall.Adna
import Bio.Iteratee
import Control.Concurrent.Async
import Data.Bits
import Data.Foldable
import Data.Ix
import Data.Maybe
import Numeric.Optimization.Algorithms.HagerZhang05
import System.Environment

import qualified Data.Vector                as V
import qualified Data.Vector.Generic        as G
import qualified Data.Vector.Unboxed        as U

import AD
import Prelude hiding ( sequence_, mapM, mapM_, concatMap, sum, minimum, foldr1 )

-- | Roughly @Maybe (Nucleotide, Nucleotide)@, encoded compactly
newtype NP = NP { unNP :: Word8 } deriving (Eq, Ord, Ix)
data Seq = Merged { unSeq :: U.Vector Word8 }
         | First  { unSeq :: U.Vector Word8 }
         | Second { unSeq :: U.Vector Word8 }

instance Show NP where
    show (NP w)
        | w  ==  16 = "NN"
        | w   >  16 = "XX"
        | otherwise = [ "ACGT" !! fromIntegral (w `shiftR` 2)
                      , "ACGT" !! fromIntegral (w .&. 3) ]


sigmoid2, isigmoid2 :: (Num a, Fractional a, Floating a) => a -> a
sigmoid2 x = y*y where y = (exp x - 1) / (exp x + 1)
isigmoid2 y = log $ (1 + sqrt y) / (1 - sqrt y)

{-# INLINE lk_fun1 #-}
lk_fun1 :: (Num a, Show a, Fractional a, Floating a, Memorable a) => Int -> [a] -> V.Vector Seq -> a
lk_fun1 lmax parms = case length parms of
    1 -> V.foldl' (\a b -> a - log (lk tab00 tab00 tab00 b)) 0 . guardV           -- undamaged case
      where
        !tab00 = fromListN (rangeSize my_bounds) [ l_epq p_subst 0 0 x
                                                 | (_,_,x) <- range my_bounds ]

    4 -> V.foldl' (\a b -> a - log (lk tabDS tabDS1 tabDS1 b)) 0 . guardV           -- double strand case
      where
        !tabDS = fromListN (rangeSize my_bounds) [ l_epq p_subst p_d p_e x
                                                 | (l,i,x) <- range my_bounds
                                                 , let p_d = mu $ lambda ^^ (1+i)
                                                 , let p_e = mu $ lambda ^^ (l-i) ]

        !tabDS1 = fromListN (rangeSize my_bounds) [ l_epq p_subst p_d 0 x
                                                  | (_,i,x) <- range my_bounds
                                                  , let p_d = mu $ lambda ^^ (1+i) ]

    5 -> V.foldl' (\a b -> a - log (lk tabSS tabSS1 tabSS2 b)) 0 . guardV           -- single strand case
      where
        !tabSS = fromListN (rangeSize my_bounds) [ l_epq p_subst p_d 0 x
                                                 | (l,i,x) <- range my_bounds
                                                 , let lam5 = lambda ^^ (1+i) ; lam3 = kappa ^^ (l-i)
                                                 , let p_d = mu $ lam3 + lam5 - lam3 * lam5 ]

        !tabSS1 = fromListN (rangeSize my_bounds) [ l_epq p_subst p_d 0 x
                                                  | (_,i,x) <- range my_bounds
                                                  , let p_d = mu $ lambda ^^ (1+i) ]

        !tabSS2 = fromListN (rangeSize my_bounds) [ l_epq p_subst 0 p_d x
                                                  | (_,i,x) <- range my_bounds
                                                  , let p_d = mu $ lambda ^^ (1+i) ]

    _ -> error "Not supposed to happen:  unexpected number of model parameters."
  where
    ~(l_subst : ~(l_sigma : ~(l_delta : ~(l_lam : ~(l_kap : _))))) = parms

    p_subst = 0.33333 * sigmoid2 l_subst
    sigma   = sigmoid2 l_sigma
    delta   = sigmoid2 l_delta
    lambda  = sigmoid2 l_lam
    kappa   = sigmoid2 l_kap

    guardV = V.filter (\u -> U.length (unSeq u) >= lmin && U.length (unSeq u) <= lmax)

    -- Likelihood given precomputed damage table.  We compute the giant
    -- table ahead of time, which maps length, index and base pair to a
    -- likelihood.
    lk tab_m     _     _ (Merged b) = U.ifoldl' (\a i np -> a * tab_m `bang` index' my_bounds (U.length b, i, NP np)) 1 b
    lk     _ tab_f     _ (First  b) = U.ifoldl' (\a i np -> a * tab_f `bang` index' my_bounds (U.length b, i, NP np)) 1 b
    lk     _     _ tab_s (Second b) = U.ifoldl' (\a i np -> a * tab_s `bang` index' my_bounds (U.length b, i, NP np)) 1 b

    index' bnds x | inRange bnds x = index bnds x
                  | otherwise = error $ "Huh? " ++ show x ++ " \\nin " ++ show bnds

    my_bounds = ((lmin,0,NP 0),(lmax,lmax,NP 16))
    mu p = sigma * p + delta * (1-p)


-- Likelihood for a certain pair of bases given error rate, C-T-rate
-- and G-A rate.
l_epq :: (Num a, Fractional a, Floating a) => a -> a -> a -> NP -> a
l_epq e p q (NP x) = case x of {
     0 -> s         ;  1 -> e         ;  2 -> e         ;  3 -> e         ;
     4 -> e         ;  5 -> s-p+4*e*p ;  6 -> e         ;  7 -> e+p-4*e*p ;
     8 -> e+q-4*e*q ;  9 -> e         ; 10 -> s-q+4*e*q ; 11 -> e         ;
    12 -> e         ; 13 -> e         ; 14 -> e         ; 15 -> s         ;
     _ -> 1 } where s = 1 - 3 * e


lkfun :: Int -> V.Vector Seq -> U.Vector Double -> Double
lkfun lmax brs parms = lk_fun1 lmax (U.toList parms) brs

combofn :: Int -> V.Vector Seq -> U.Vector Double -> (Double, U.Vector Double)
combofn lmax brs parms = (x,g)
  where D x g = lk_fun1 lmax (paramVector $ U.toList parms) brs

params :: Parameters
params = defaultParameters { printFinal = False, verbose = Quiet, maxItersFac = 20 }

lmin :: Int
lmin = 25

main :: IO ()
main = do
    [fp] <- getArgs
    brs <- subsampleBam fp >=> run $ \_ ->
           joinI $ filterStream (\b -> not (isUnmapped (unpackBam b)) && G.length (b_seq (unpackBam b)) >= lmin) $
           joinI $ takeStream 100000 $
           joinI $ mapStream pack_record $
           joinI $ filterStream (\u -> U.length (U.filter (<16) (unSeq u)) * 10 >= 9 * U.length (unSeq u)) $
           stream2vectorN 30000

    let lmax = V.maximum $ V.map (U.length . unSeq) brs
        v0 = crude_estimate brs
        opt v = optimize params 0.0001 v
                         (VFunction $ lkfun lmax brs)
                         (VGradient $ snd . combofn lmax brs)
                         (Just . VCombined $ combofn lmax brs)

    results <- mapConcurrently opt [ v0, U.take 4 v0, U.take 1 v0 ]

    let mlk = minimum [ finalValue st | (_,_,st) <- results ]
        tot = sum [ exp $ mlk - finalValue st | (_,_,st) <- results ]
        p l = exp (mlk - l) / tot

        [ (p_ss, [ _, ssd_sigma_, ssd_delta_, ssd_lambda, ssd_kappa ]),
          (p_ds, [ _, dsd_sigma_, dsd_delta_, dsd_lambda ]),
          (_   , [ _ ]) ] = [ (p (finalValue st), map sigmoid2 $ G.toList xs) | (xs,_,st) <- results ]

        ssd_sigma = p_ss * ssd_sigma_
        ssd_delta = p_ss * ssd_delta_
        dsd_sigma = p_ds * dsd_sigma_
        dsd_delta = p_ds * dsd_delta_

    print DP{..}

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
pack_record br = if isReversed b then k (revcom u1) else k u1
  where
    b@BamRec{ b_cigar = Cigar cigar, .. } = unpackBam br

    k | isMerged     b = Merged
      | isTrimmed    b = Merged
      | isSecondMate b = Second
      | otherwise      = First

    revcom = U.reverse . U.map (\x -> if x > 15 then x else xor x 15)
    u1 = U.fromList . map unNP $ go cigar (G.toList b_seq) (fromMaybe [] $ getMd b)

    go :: [(CigOp,Int)] -> [Nucleotides] -> [MdOp] -> [NP]

    go ((_,0):cs)   ns mds  = go cs ns mds
    go cs ns (MdNum  0:mds) = go cs ns mds
    go cs ns (MdDel []:mds) = go cs ns mds
    go  _ []              _ = []
    go []  _              _ = []

    go ((Mat,nm):cs) (n:ns) (MdNum mm:mds) = mk_pair n n  : go ((Mat,nm-1):cs) ns (MdNum (mm-1):mds)
    go ((Mat,nm):cs) (n:ns) (MdRep n':mds) = mk_pair n n' : go ((Mat,nm-1):cs) ns               mds
    go ((Mat,nm):cs)    ns  (MdDel _ :mds) =                go ((Mat, nm ):cs) ns               mds

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
    mk_pair' u (Ns b) = case b of 1 -> NP $ u .|. 0
                                  2 -> NP $ u .|. 4
                                  4 -> NP $ u .|. 8
                                  8 -> NP $ u .|. 12
                                  _ -> esc


infix 7 /%/
(/%/) :: Integral a => a -> a -> Double
0 /%/ 0 = 0
a /%/ b = fromIntegral a / fromIntegral b

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
-- undamaged fraction:  see below                                      √
--
-- Contaminant fraction:  let f5 (f3, f1) be the fraction of reads
-- showing damage at the 5' end (3' end, both ends).  Let a (b) be
-- the probability of an endogenous reads to show damage at the 5'
-- end (3' end).  Let e be the fraction of endogenous reads.  Then
-- we have:
--
-- f5 = e * a
-- f3 = e * b
-- f1 = e * a * b
--
-- f5 * f3 / f1 = e
--
-- Straight forward and easy to understand, but in practice, this method
-- produces ridiculous overestimates, ridiculous underestimates,
-- negative contamination rates, and general grief.  It's actually
-- better to start from a constant number.


crude_estimate :: V.Vector Seq -> U.Vector Double
crude_estimate seqs0 = U.fromList [ l_subst, l_sigma, l_delta, l_lam, l_kap ]
  where
    seqs = V.filter ((>= 10) . U.length) $ V.map unSeq seqs0

    total_equals = V.sum (V.map (U.length . U.filter      isNotSubst) seqs)
    total_substs = V.sum (V.map (U.length . U.filter isOrdinarySubst) seqs) * 6 `div` 5
    l_subst = isigmoid2 $ max 0.001 $ total_substs /%/ (total_equals + total_substs)

    c_to_t, g_to_a, c_to_c :: Word8
    c_to_t = 7
    g_to_a = 8
    c_to_c = 5

    isNotSubst x = x < 16 && x `shiftR` 2 == x .&. 3
    isOrdinarySubst x = x < 16 && x `shiftR` 2 /= x .&. 3 &&
                        x /= c_to_t && x /= g_to_a

    ct_at_alpha = V.length $ V.filter (\v -> v U.! 0 == c_to_t && dmg_omega v) seqs
    cc_at_alpha = V.length $ V.filter (\v -> v U.! 0 == c_to_c && dmg_omega v) seqs
    ct_at_beta  = V.length $ V.filter (\v -> v U.! 1 == c_to_t && dmg_omega v) seqs
    cc_at_beta  = V.length $ V.filter (\v -> v U.! 1 == c_to_c && dmg_omega v) seqs

    dmg_omega v = v U.! (l-1) == c_to_t || v U.! (l-1) == g_to_a
               || v U.! (l-2) == c_to_t || v U.! (l-2) == g_to_a
               || v U.! (l-3) == c_to_t || v U.! (l-3) == g_to_a
        where l = U.length v

    l_lam = isigmoid2 lambda
    lambda = min 0.9 $ max 0.1 $
                (ct_at_beta * (cc_at_alpha + ct_at_alpha)) /%/
                ((cc_at_beta + ct_at_beta) * ct_at_alpha)

    ct_at_omega = V.length $ V.filter (\v -> v U.! (U.length v -1) == c_to_t && dmg_alpha v) seqs
    cc_at_omega = V.length $ V.filter (\v -> v U.! (U.length v -1) == c_to_c && dmg_alpha v) seqs
    ct_at_psi   = V.length $ V.filter (\v -> v U.! (U.length v -2) == c_to_t && dmg_alpha v) seqs
    cc_at_psi   = V.length $ V.filter (\v -> v U.! (U.length v -2) == c_to_c && dmg_alpha v) seqs

    dmg_alpha v = v U.! 0 == c_to_t || v U.! 1 == c_to_t || v U.! 2 == c_to_t

    l_kap = isigmoid2 $ min 0.9 $ max 0.1 $
                (ct_at_psi * (cc_at_omega+ct_at_omega)) /%/
                ((cc_at_psi+ct_at_psi) * ct_at_omega)

    total_inner_CCs = V.sum $ V.map (U.length . U.filter (== c_to_c) . takeInner) seqs
    total_inner_CTs = V.sum $ V.map (U.length . U.filter (== c_to_t) . takeInner) seqs
    takeInner v = U.slice 5 (U.length v - 10) v

    delta = (total_inner_CTs /%/ (total_inner_CTs+total_inner_CCs))
    raw_rate = ct_at_alpha /%/ (ct_at_alpha + cc_at_alpha)

    -- clamping is necessary if f_endo ends up wrong
    l_delta = isigmoid2 $ min 0.99 delta
    l_sigma = isigmoid2 . min 0.99 $ raw_rate / lambda


class Memorable a where
    type Memo a :: *

    fromListN :: Int -> [a] -> Memo a
    bang :: Memo a -> Int -> a

instance Memorable Double where
    type Memo Double = U.Vector Double

    fromListN = U.fromListN
    bang = (U.!)

instance Memorable AD where
    type Memo AD = (Int, U.Vector Double)

    fromListN n xs@(D _ v:_) = (1+d, U.fromListN (n * (1+d)) $ concatMap unpack xs)
      where
        !d = U.length v
        unpack (C a)    = a : replicate d 0
        unpack (D a da) = a : U.toList da

    bang (d, v) i = D (v U.! (d*i+0)) (U.slice (d*i+1) (d-1) v)
