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
--   - needs better output
--   - needs support for multiple input files
--   - needs to deal with long (unmerged) reads (by ignoring them?)

import Bio.Adna              hiding ( bang )
import Bio.Bam.Header
import Bio.Bam.Index
import Bio.Bam.Rec
import Bio.Base
import Bio.Genocall.Metadata
import Bio.Iteratee
import Bio.Util.AD
import Bio.Util.AD2
import Bio.Util.Numeric hiding ( sigmoid2, isigmoid2 )
import Control.Applicative
import Control.Concurrent.Async
import Control.Monad                ( unless )
import Data.Bits
import Data.Foldable
import Data.Ix
import Data.Maybe
import Data.String                  ( fromString )
import Data.Text                    ( unpack )
import System.Console.GetOpt
import System.Environment
import System.Exit
import System.FilePath
import System.IO                    ( hPutStrLn )

import qualified Data.HashMap.Strict        as M
import qualified Data.Vector                as V
import qualified Data.Vector.Generic        as G
import qualified Data.Vector.Unboxed        as U

import Prelude hiding ( sequence_, mapM, mapM_, concatMap, sum, minimum, foldr1, foldl )

-- | Full likelihood function.  We apply the HKY substitution model
-- (anything simpler doesn't seem to cut it), and we will estimate a
-- position specific substitution rate for C->T and G->A.  We
-- arbitrarily estimate (k+1) such values (say 12) for each end, and another
-- one for the interior.
--
-- This leaves us with the following parameters:  pi_a, pi_c, pi_g for
-- the base composition; mu and nu for the rates of transitions and
-- transversions, respectively; alpha_d and beta_d for the C->T and G->A
-- substitution rates in the doubles stranded interior; alpha_{0..k} and
-- beta_{0..k} for the substitution rates at positions from the 5' end;
-- alpha'_{0..k} and beta'_{0..k} for the substitution rates at
-- positions from the 3' end.
--
-- The remaining input is two sets of 'SubstitutionStats', one each for
-- 5' and 3' ends.
{-# INLINE lk_fun #-}
lk_fun :: (Num a, Show a, Fractional a, Floating a, Memorable a)
        => SubstitutionStats -> SubstitutionStats -> [a] -> a
lk_fun stats5 stats3 (pi_a : pi_c : pi_g : mu0 : nu0 : alpha : beta : spars) =
    [
    -- pick vectors of counts for each of the possible substitutions
    | fromNuc <- [nucA,nucC,nucG,nucT] | row <- smat
    ???
  where
    pp = recip $ exp pi_a + exp pi_c + exp pi_g + 1
    pA = exp pi_a * pp
    pC = exp pi_c * pp
    pG = exp pi_g * pp
    pT =            PP

    mu = recip $ 1 + exp (-mu0)
    nu = recip $ 1 + exp (-nu0)

    smat alpha0 beta0 =
        [ [ 1 - nu*pC - mu*pG - nu*pT + mu*pG * beta,
                            nu*pC * alpham1,
                                            mu*pG * betam1,
                                                            nu*pC * alpha + nu*pT ]
        , [ nu*pA + nu*pG * beta,
                            (1 - nu*pA - nu*pG - mu*pT) * alpham1,
                                            nu*pG * betam1,
                                                            (1 - nu*pA - nu*pG - mu*pT) * alpha + mu * pT ]
        , [ (1 - mu*pA - nu*pC - nu_pT) * beta + mu*pA,
                            nu*pC * alpham1,
                                            (1 - mu*pA - nu*pC - nu*pT) * betam1,
                                                            nu*pC * alpha + nu*pT ],
        , [ nu*pA + nu*pG * beta,
                            mu*pC * alpham1,
                                            mu*pG * betam1,
                                                            1 - nu*pA - mu*pC - nu*pG + mu*pC * alpha ] ]
      where
        alpham1 = recip $ 1 + exp alpha0
        alpha   = recip $ 1 + exp (-alpha0)
        betam1  = recip $ 1 + exp beta0
        beta    = recip $ 1 + exp (-beta0)






case length parms of
    1 -> V.foldl' (\a b -> a - log (lk tab00 tab00 tab00 b)) 0 . guardV             -- undamaged case
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
    lk tab_m     _     _ (Merged  b) = U.ifoldl' (\a i np -> a * tab_m `bang` index' my_bounds (U.length b, i, NP np)) 1 b
    lk     _ tab_f     _ (Mate1st b) = U.ifoldl' (\a i np -> a * tab_f `bang` index' my_bounds (U.length b, i, NP np)) 1 b
    lk     _     _ tab_s (Mate2nd b) = U.ifoldl' (\a i np -> a * tab_s `bang` index' my_bounds (U.length b, i, NP np)) 1 b

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


-- Likelihood for a certain pair of bases given error rate, C-T-rate
-- and G-A rate.
l_abp :: (Num a, Fractional a, Floating a) => a -> a -> a -> NP -> a
l_abp e p q (NP x) = case x of {
     0 -> s         ;  1 -> e         ;  2 -> e         ;  3 -> e         ;
     4 -> e         ;  5 -> s-p+4*e*p ;  6 -> e         ;  7 -> e+p-4*e*p ;
     8 -> e+q-4*e*q ;  9 -> e         ; 10 -> s-q+4*e*q ; 11 -> e         ;
    12 -> e         ; 13 -> e         ; 14 -> e         ; 15 -> s         ;
     _ -> 1 } where s = 1 - 3 * e


lkfun :: Int -> Int -> V.Vector Seq -> U.Vector Double -> Double
lkfun lmin lmax brs parms = lk_fun1 lmin lmax (U.toList parms) brs

lkfun' :: Int -> Int -> V.Vector Seq -> [Double] -> AD
lkfun' lmin lmax brs parms = lk_fun1 lmin lmax (paramVector parms) brs

lkfun'' :: Int -> Int -> V.Vector Seq -> [Double] -> AD2
lkfun'' lmin lmax brs parms = lk_fun1 lmin lmax (paramVector2 parms) brs

combofn :: Int -> Int -> V.Vector Seq -> U.Vector Double -> (Double, U.Vector Double)
combofn lmin lmax brs parms = (x,g)
  where D x g = lk_fun1 lmin lmax (paramVector $ U.toList parms) brs


data Conf = Conf {
    conf_lmin :: Int,
    conf_metadata :: FilePath,
    conf_leehom :: Bool,
    conf_report :: String -> IO (),
    conf_params :: Parameters }

defaultConf :: Conf
defaultConf = Conf 25 (error "no config file specified") False (\_ -> return ()) quietParameters

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "m"  ["min-length"]   (ReqArg set_lmin  "LEN") "Set minimum length to LEN (25)",
    Option "c"  ["config"]       (ReqArg set_conf "FILE") "Configuiration is stored in FILE",
    Option [ ]  ["leehom"]       (NoArg       set_leehom) "Flags are broken because of leeHom",
    Option "v"  ["verbose"]      (NoArg      set_verbose) "Print progress reports",
    Option "h?" ["help","usage"] (NoArg       disp_usage) "Print this message and exit" ]
  where
    set_lmin   a c = readIO a >>= \l -> return $ c { conf_lmin     = l }
    set_conf   f c =                    return $ c { conf_metadata = f }
    set_leehom f c =                    return $ c { conf_leehom   = True }
    set_verbose  c = return $ c { conf_report   = hPutStrLn stderr, conf_params = debugParameters }

    disp_usage  _ = do pn <- getProgName
                       let blah = "Usage: " ++ pn ++ " [OPTION...] [LIBRARY-NAME...]"
                       putStrLn $ usageInfo blah options
                       exitSuccess

main :: IO ()
main = do
    (opts, lnames, errors) <- getOpt Permute options <$> getArgs
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    conf <- foldl (>>=) (return defaultConf) opts
    mapM_ (main' conf) lnames

main' :: Conf -> String -> IO ()
main' Conf{..} lname = do
    [Library _ fs _] <- return . filter ((fromString lname ==) . library_name) . concatMap sample_libraries . M.elems
                        =<< readMetadata conf_metadata

    -- XXX  meh.  subsampling from multiple files is not yet supported :(
    -- We subsample up to a million reads; that should be plenty, and it
    -- should also be plenty fast now.  All we need is substitution
    -- statistics, and we'll collect them, say, 50 positions into the
    -- read.  That should be plenty.

    DmgStats{substs5,substs3} <-
            subsampleBam (takeDirectory conf_metadata </> unpack (head fs)) >=> run $ \_ ->
            joinI $ filterStream (\b -> not (isUnmapped (unpackBam b)) && G.length (b_seq (unpackBam b)) >= conf_lmin) $
            joinI $ takeStream 1000000 $
            joinI $ damagePatternsIterMD 50 conf_leehom $ skipToEof

    let lmax = V.maximum $ V.map (U.length . unSeq) brs
        v0 = crude_estimate brs
        opt v = optimize conf_params 0.0001 v
                         (VFunction $ lkfun conf_lmin lmax brs)
                         (VGradient $ snd . combofn conf_lmin lmax brs)
                         (Just . VCombined $ combofn conf_lmin lmax brs)

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

    putStrLn $ "p_{ss} = " ++ show p_ss ++ ", p_{ds} = " ++ show p_ds
    putStrLn $ show DP{..}
    updateMetadata (store_dp lname DP{..}) conf_metadata

    -- Trying to get confidence intervals.  Right now, just get the
    -- gradient and Hessian at the ML point.  Gradient should be nearly
    -- zero, Hessian should be symmetric and positive definite.
    -- (Remember, we minimized.)
    mapM_ print [ (r,s) | (_,r,s) <- results ]
    putStrLn ""
    mapM_ print [ lkfun' conf_lmin lmax brs (G.toList xs) | (xs,_,_) <- results ]
    putStrLn ""
    mapM_ print [ lkfun'' conf_lmin lmax brs (G.toList xs) | (xs,_,_) <- results ]

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

    fromListN _    [       ] = error "unexpected: tried to memorize an empty list"
    fromListN _    (C _  :_) = error "unexpected: tried to memorize a value without derivatives"
    fromListN n xs@(D _ v:_) = (1+d, U.fromListN (n * (1+d)) $ concatMap unp xs)
      where
        !d = U.length v
        unp (C a)    = a : replicate d 0
        unp (D a da) = a : U.toList da

    bang (d, v) i = D (v U.! (d*i+0)) (U.slice (d*i+1) (d-1) v)

instance Memorable AD2 where
    type Memo AD2 = (Int, U.Vector Double)

    fromListN _    [            ] = error "unexpected: tried to memorize an empty list"
    fromListN _    (C2 _     : _) = error "unexpected: tried to memorize a value without derivatives"
    fromListN n xs@(D2 _ v _ : _) = (d, U.fromListN (n * (1+d+d*d)) $ concatMap unp xs)
      where
        !d = U.length v
        unp (C2 a)        = a : replicate (d+d*d) 0
        unp (D2 a da dda) = a : U.toList da ++ U.toList dda

    bang (d, v) i = D2 (v U.! (stride*i))
                       (U.slice (stride*i+1) d v)
                       (U.slice (stride*i+1+d) (d*d) v)
      where
        stride = 1 + d + d*d


store_dp :: String -> DamageParameters Double -> Metadata -> Metadata
store_dp lname dp = M.map go1
  where
    go1 (Sample  ls af bf ts dv) = Sample (map go2 ls) af bf ts dv
    go2 (Library      nm fs dmg)
        | nm == fromString lname = Library nm fs (Just dp)
        | otherwise              = Library nm fs dmg



sigmoid2 x = recip $ 1 + exp (1-x)
isigmoid2 y = 1 - log ((1-y) / y)
