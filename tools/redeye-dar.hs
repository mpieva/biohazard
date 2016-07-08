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
import Bio.Bam.Index
import Bio.Bam.Rec
import Bio.Base
import Bio.Genocall.Metadata
import Bio.Iteratee
import Bio.Util.AD
import Bio.Util.AD2
import Control.Applicative
import Control.Monad                ( unless )
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
import qualified Data.Vector.Generic        as V
import qualified Data.Vector.Unboxed        as U

import Prelude hiding ( sequence_, mapM, mapM_, concatMap, sum, minimum, foldr1, foldl )

-- | Full likelihood function.  We apply the HKY substitution model
-- (anything simpler doesn't seem to cut it), and we will estimate a
-- position specific substitution rate for C->T and G->A.  We
-- arbitrarily estimate (k+1) such values (more than 8 doesn't look very
-- useful) for each end, and another one for the interior.
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
--
-- Notes:
--
-- - The mutation rates could add up to more than one.  Then the
--   likelihood becomes NaN.  We might try to reformulate to avoid this.
--
-- - We could compute the full likelihood function with gradient and
--   Hessian once at the end, just to see if the fit is good and how
--   badly the parameters are confounded.
--
-- - Alternatively, we could parameterize the deamination rates using an
--   exponential decay.  The old formulation had three models:  one
--   fixes deamination to zero, one has C->T at one end and G->A at the
--   other with one length parameter, the third has only C->T with two
--   length parameters.  This works out, respectively, to
--
--      alpha5_i = alpha5_0 * lambda ^^ (1+i)
--      beta5_i  = 0
--      alpha3_i = 0
--      beta3_i  = beta3_0 * lambda ^^ (1+i)
--
--   and
--
--      alpha5_i = alpha5_0 * lambda ^^ (1+i)
--      beta5_i  = 0
--      alpha3_i = alpha3_0 * kappa ^^ (1+i)
--      beta3_i  = 0-

-- Argh, dumb mistake.  Likelihood of any substitution has to be
-- multiplied with base frequency.
smat :: (Num a, Floating a, IsDouble a) => Double -> a -> a -> a -> a -> [a]
smat gc mu nu a b =
       [ pA * (1 - nu*pC - mu*pG - nu*pT + mu*pG * tr (-b))
                   , pA * (nu*pC * tr a)
                               , pA * (mu*pG * tr b)
                                           , pA * (nu*pC * tr (-a) + nu*pT)
       , pC * (nu*pA + nu*pG * tr (-b))
                    , pC * ((1 - nu*pA - nu*pG - mu*pT) * tr a)
                               , pC * (nu*pG * tr b)
                                           , pC * ((1 - nu*pA - nu*pG - mu*pT) * tr (-a) + mu * pT)
       , pG * ((1 - mu*pA - nu*pC - nu*pT) * tr (-b) + mu*pA)
                    , pG * (nu*pC * tr a)
                               , pG * ((1 - mu*pA - nu*pC - nu*pT) * tr b)
                                           , pG * (nu*pC * tr (-a) + nu*pT)
       , pT * (nu*pA + nu*pG * tr (-b))
                    , pT * (mu*pC * tr a)
                               , pT * (mu*pG * tr b)
                                           , pT * (1 - nu*pA - mu*pC - nu*pG + mu*pC * tr (-a)) ]
  where
    tr x = recip $ 1 + exp x

    pA = fromDouble $ 0.5 - gc * 0.5
    pC = fromDouble $ gc * 0.5
    pG = fromDouble $ gc * 0.5
    pT = fromDouble $ 0.5 - gc * 0.5


-- | Estimates GC fraction.  Base composition is stable with respect to
-- our substitution model, so we should be able to estimate it from the
-- reference sequence.  Practically, since G and C have to be balanced,
-- we just estimate the GC fraction.

est_gc_frac :: DmgStats a -> Double
est_gc_frac DmgStats{..} = fromIntegral gc_counts / fromIntegral (at_counts+gc_counts)
  where
    get n = maybe 0 U.sum . lookup (Just n)
    at_counts = sum [ get n v | n <- [nucA,nucT], v <- [basecompo5, basecompo3] ]
    gc_counts = sum [ get n v | n <- [nucG,nucC], v <- [basecompo5, basecompo3] ]

-- | Estimates mutation rates for transitions and transversion.  We
-- write down the likelihood function with just two general deamination
-- rates.  Likewise, our substitution statistics collapse to one single
-- matrix.  This is not perfect, but hopefully good enough to get a
-- quick estimate of the two mutation rates.  (Numeric optimization,
-- sorry it's in IO.)
est_mut_rate :: Parameters -> Double -> DmgStats a -> IO (Double, Double)
est_mut_rate params gc_frac DmgStats{..} = do
    (xs,_,_) <- minimize params 0.001 lk_fun v0
    return (xs V.! 0, xs V.! 1)
  where
    counts = [ c5 + c3
             | subst <- range (nucA :-> nucA, nucT :-> nucT)
             , let c5 = maybe 0 U.sum $ lookup subst substs5
             , let c3 = maybe 0 U.sum $ lookup subst substs3 ]

    v0 = U.fromList [-4.5, -5.5, -3, -6.5] -- guesstimates

    lk_fun [mu0,nu0,alpha,beta] =
        let !mu = recip $ 1 + exp (-mu0)
            !nu = recip $ 1 + exp (-nu0)
        in negate . sum $ zipWith (\num lk -> fromIntegral num * log lk)
                                  counts (smat gc_frac mu nu alpha beta)

-- | Estimate deamination rates.  We take the base composition (GC
-- fraction) and the two mutation rates (ti,tv) as constant, then
-- estimate two deamination rates (C->T and G->A) for a given number of
-- positions, and then the remainder.  (They are now decoupled, so we
-- could as well estimate the parameters (just two!) using Newton
-- iteration.)
est_deamination :: Parameters -> Int -> Double -> (Double,Double) -> DmgStats a
                -> IO (U.Vector Double, Double, Double, U.Vector Double)
est_deamination parameters l_over gc_frac (mu0,nu0) DmgStats{..} = do
    ab_left <- V.concat <$> sequence
               [ fst3 <$> minimize parameters 0.001 (lk_fun counts) v0
               | i <- [0 .. l_over-1]
               , let counts = [ maybe 0 (U.! i) $ lookup subst substs5
                              | subst <- range (nucA :-> nucA, nucT :-> nucT) ] ]

    ab_mid <- let counts = [ maybe 0 (U.sum . U.drop l_over)             (lookup subst substs5)
                           + maybe 0 (U.sum . U.drop l_over . U.reverse) (lookup subst substs3)
                           | subst <- range (nucA :-> nucA, nucT :-> nucT) ]
              in fst3 <$> minimize parameters 0.001 (lk_fun counts) v0

    ab_right <- V.concat <$> sequence
                [ fst3 <$> minimize parameters 0.001 (lk_fun counts) v0
                | i <- [0 .. l_over-1]
                , let counts = [ maybe 0 ((U.! i) . U.reverse) $ lookup subst substs3
                               | subst <- range (nucA :-> nucA, nucT :-> nucT) ] ]

    return (V.convert ab_left, ab_mid V.! 0, ab_mid V.! 1, V.convert ab_right)

  where
    fst3 (x,_,_) = x
    !mu = fromDouble . recip $ 1 + exp (-mu0)
    !nu = fromDouble . recip $ 1 + exp (-nu0)

    v0 = U.fromList [-3, -6.5] -- guesstimates

    lk_fun :: [Int] -> [AD] -> AD
    lk_fun counts [alpha,beta] =
        negate . sum $ zipWith (\num lk -> fromIntegral num * log lk)
                               counts (smat gc_frac mu nu alpha beta)


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
    set_leehom   c =                    return $ c { conf_leehom   = True }
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

    dmg_stats <- subsampleBam (takeDirectory conf_metadata </> unpack (head fs)) >=> run $ \_ ->
                 joinI $ filterStream (\b -> not (isUnmapped (unpackBam b)) && V.length (b_seq (unpackBam b)) >= conf_lmin) $
                 joinI $ takeStream 100000 $
                 damagePatternsIterMD 50 conf_leehom skipToEof

    let gc_frac = est_gc_frac dmg_stats
    print gc_frac

    (mu,nu) <- est_mut_rate conf_params gc_frac dmg_stats
    print (mu,nu)

    (ab_left, alpha_mid, beta_mid, ab_right) <- est_deamination conf_params 12 gc_frac (mu,nu) dmg_stats
    print ab_left
    print (alpha_mid, beta_mid)
    print ab_right

    {-
    let ct i = log 0.5 + log 0.6 * i
        ga i = log 0.05 + log 0.6 * i
        omax = 8

    let v0 = 0 : -0.4 : -0.4 : -4.5 : -6 : ct omax : ga omax : concat
             [ [ ct i, ga i, ct i, ga i ] | i <- [0 .. omax-1] ]
    results <- minimize conf_params 0.001 (lk_fun substs5 substs3) (U.fromList v0)
                         -- (VFunction $ lkfun substs5 substs3)
                         -- (VGradient $ snd . combofn substs5 substs3)
                         -- (Just . VCombined $ combofn substs5 substs3)
    print results
    -}
    --
    -- results <- mapConcurrently opt [ v0, U.take 4 v0, U.take 1 v0 ]

    -- let mlk = minimum [ finalValue st | (_,_,st) <- results ]
        -- tot = sum [ exp $ mlk - finalValue st | (_,_,st) <- results ]
        -- p l = exp (mlk - l) / tot

        {-[ (p_ss, [ _, ssd_sigma_, ssd_delta_, ssd_lambda, ssd_kappa ]),
          (p_ds, [ _, dsd_sigma_, dsd_delta_, dsd_lambda ]),
          (_   , [ _ ]) ] = [ (p (finalValue st), map sigmoid2 $ V.toList xs) | (xs,_,st) <- results ]

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
    mapM_ print [ lkfun' conf_lmin lmax brs (V.toList xs) | (xs,_,_) <- results ]
    putStrLn ""
    mapM_ print [ lkfun'' conf_lmin lmax brs (V.toList xs) | (xs,_,_) <- results ] -}


store_dp :: String -> DamageParameters Double -> Metadata -> Metadata
store_dp lname dp = M.map go1
  where
    go1 (Sample  ls af bf ts dv) = Sample (map go2 ls) af bf ts dv
    go2 (Library      nm fs dmg)
        | nm == fromString lname = Library nm fs (OldDamage dp)
        | otherwise              = Library nm fs dmg
