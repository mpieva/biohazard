{-# LANGUAGE RecordWildCards #-}
module Bio.Genocall.Damage where

import Bio.Adna
import Bio.Bam
import Bio.Prelude
import Bio.Util.AD

import qualified Data.Vector.Generic as V
import qualified Data.Vector.Unboxed as U

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

estimateDamage :: Parameters -> DmgStats a -> IO (NewDamageParameters U.Vector Double)
estimateDamage conf_params dmg = do
    let dp_gc_frac = est_gc_frac dmg
    (dp_mu,dp_nu) <- est_mut_rate conf_params dp_gc_frac dmg

    (ab_left, dp_alpha, dp_beta, ab_right) <- est_deamination conf_params 12 dp_gc_frac (dp_mu,dp_nu) dmg
    let dp_alpha5 = U.ifilter (const .  odd) ab_left
        dp_beta5  = U.ifilter (const . even) ab_left
        dp_alpha3 = U.ifilter (const .  odd) ab_right
        dp_beta3  = U.ifilter (const . even) ab_right
    return NDP{..}

estimateDamageFromFiles :: Int -> Parameters -> [String] -> IO (NewDamageParameters U.Vector Double)
estimateDamageFromFiles lmin params fs = do
    foldM (\acc f -> decodeAnyBamFile f >=> run >=> return . mappend acc $ \hdr ->
                     filterStream (\b ->
                            not (isUnmapped (unpackBam b)) &&
                            V.length (b_seq (unpackBam b)) >= lmin)     =$
                     takeStream 100000                                  =$
                     damagePatternsIterMD 50 (leehom hdr) skipToEof)
          mempty fs
    >>= estimateDamage params
  where
    leehom :: BamMeta -> Bool
    leehom meta = not $ null [ () | ("PG",line) <- meta_other_shit meta
                                  , ("PN","mergeTrimReadsBAM") <- line ]
