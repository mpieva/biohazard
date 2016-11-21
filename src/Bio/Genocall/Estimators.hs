{-# LANGUAGE DeriveGeneric #-}
{-# OPTIONS_GHC -O0 #-}
module Bio.Genocall.Estimators (
        estimateDamageFromFiles,
        tabulateSingle,
        estimateSingle,
        DivEst(..),
        DivTable(..),
        good_regions
    ) where

-- XXX  Many optimizations fit only one parameter.  Newton-Iteration
-- should be more efficient than the generic CG method.
--
-- For Newton-Iteration, the update is 'x := x + dt' where
-- 'dt = -H(x)^{-1} \Nabla f(x)' or 'H(x) dt = -\Nabla f(x)'

import Bio.Adna
import Bio.Bam
import Bio.Bam.Pileup                ( p_snp_pile )
import Bio.Genocall                  ( Snp_GLs(..), Calls )
import Bio.Prelude
import Bio.Util.AD
import Bio.Util.AD2
import Bio.Util.Numeric              ( log1pexp )
import Bio.Util.Pretty
import Data.Binary
import Data.Vector.Binary            ()

import qualified Data.Vector.Generic as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as M

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
-- * The mutation rates could add up to more than one.  Then the
--   likelihood becomes NaN.  We might try to reformulate to avoid this.
--
-- * We could compute the full likelihood function with gradient and
--   Hessian once at the end, just to see if the fit is good and how
--   badly the parameters are confounded.
--
-- * Alternatively, we could parameterize the deamination rates using an
--   exponential decay.  The old formulation had three models:  one
--   fixes deamination to zero, one has C->T at one end and G->A at the
--   other with one length parameter, the third has only C->T with two
--   length parameters.  This works out, respectively, to
-- \[
--      alpha5_i = alpha5_0 * lambda ^^ (1+i)
--      beta5_i  = 0
--      alpha3_i = 0
--      beta3_i  = beta3_0 * lambda ^^ (1+i)
-- \]
--   and
-- \[
--      alpha5_i = alpha5_0 * lambda ^^ (1+i)
--      beta5_i  = 0
--      alpha3_i = alpha3_0 * kappa ^^ (1+i)
--      beta3_i  = 0-
-- \]
--
-- * Dumb mistake to be avoided:  The likelihood of any substitution has
--   to be multiplied with the base frequency of the original base.

smat :: (Floating a, IsDouble a) => Double -> a -> a -> a -> a -> [a]
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
    getV n = maybe 0 U.sum . lookup (Just n)
    at_counts = sum [ getV n v | n <- [nucA,nucT], v <- [basecompo5, basecompo3] ]
    gc_counts = sum [ getV n v | n <- [nucG,nucC], v <- [basecompo5, basecompo3] ]

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

    lk_fun [mu0,nu0,alpha,beta] = sum $ zipWith (\num lk -> - fromIntegral num * log lk)
                                                counts (smat gc_frac mu nu alpha beta)
      where !mu = recip $ 1 + exp (-mu0)
            !nu = recip $ 1 + exp (-nu0)
    lk_fun _ = error "Inconceivable!"


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
    lk_fun _ _ = error "Inconceivable!"

estimateDamage :: Parameters -> DmgStats a -> IO (NewDamageParameters U.Vector Double)
estimateDamage conf_params dmg = do
    let dp_gc_frac = est_gc_frac dmg
    (dp_mu,dp_nu) <- est_mut_rate conf_params dp_gc_frac dmg

    (ab_left, dp_alpha, dp_beta, ab_right) <- est_deamination conf_params 8 dp_gc_frac (dp_mu,dp_nu) dmg
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
                     addFragType hdr                                    =$
                     damagePatternsIterMD 50 skipToEof)
          mempty fs
    >>= estimateDamage params

data DivTable = DivTable !Double !(U.Vector Int) deriving (Show, Generic)

instance Binary DivTable
instance Pretty DivTable where pretty = default_pretty
instance Parse  DivTable where parse  = default_parse

instance Monoid DivTable where
    mempty = DivTable 0 U.empty

    DivTable a u `mappend` DivTable b v = DivTable (a+b) w
      where w | U.null  u = v
              | U.null  v = u
              | otherwise = U.zipWith (+) u v

-- | Divergence estimate.  Lists contain three or four floats, these are
-- divergence, heterozygosity at W sites, heterozygosity at S sites, and
-- optionally gappiness in this order.
data DivEst = DivEst {
    point_est :: [Double],
    conf_region :: [( [Double], [Double] )]
  } deriving (Show, Generic)

instance Pretty DivEst where pretty = default_pretty
instance Parse  DivEst where parse  = default_parse

-- XXX we should estimate an indel rate, to be appended as the fourth
-- result (but that needs different tables)
estimateSingle :: DivTable -> IO (DivEst, DivEst)
estimateSingle (DivTable _llk_rr tab) = do
    (fit1, _res1, _stats1) <- minimize quietParameters 0.0001 (llk tab) (U.fromList   [0,0])
    (fit2, _res2, _stats2) <- minimize quietParameters 0.0001 (llk tab) (U.fromList [0,0,0])

    let xform v = map (\x -> recip $ 1 + exp (-x)) $ V.toList v
        !de1 = DivEst (xform fit1) (map (xform *** xform) $ confidenceIntervals (llk2 tab) (V.convert fit1))
        !de2 = DivEst (xform fit2) (map (xform *** xform) $ confidenceIntervals (llk2 tab) (V.convert fit2))

    return (de1,de2)

llk :: U.Vector Int -> [AD] -> AD
llk tab [delta,eta] = llk' tab 0 delta eta + llk' tab 6 delta eta
llk tab [delta,eta,eta2] = llk' tab 0 delta eta + llk' tab 6 delta eta2
llk _ _ = error "Wtf? (3)"

llk2 :: U.Vector Int -> [AD2] -> AD2
llk2 tab [delta,eta] = llk' tab 0 delta eta + llk' tab 6 delta eta
llk2 tab [delta,eta,eta2] = llk' tab 0 delta eta + llk' tab 6 delta eta2
llk2 _ _ = error "Wtf? (4)"

{-# INLINE llk' #-}
llk' :: (Ord a, Floating a) => U.Vector Int -> Int -> a -> a -> a
llk' tab base delta eta = block (base+0) g_RR g_RA g_AA
                        + block (base+1) g_RR g_AA g_RA
                        + block (base+2) g_RA g_RR g_AA
                        + block (base+3) g_RA g_AA g_RR
                        + block (base+4) g_AA g_RR g_RA
                        + block (base+5) g_AA g_RA g_RR
  where
    !maxD2 = U.length tab `div` 12
    !maxD  = round (sqrt (fromIntegral maxD2) :: Double)

    !g_AA = Pr delta / Pr (log1pexp delta)
    !g_RA =        1 / Pr (log1pexp delta) * Pr eta / Pr (log1pexp eta)
    !g_RR =        1 / Pr (log1pexp delta) *      1 / Pr (log1pexp eta)

    block ix g1 g2 g3 = U.ifoldl' step 0 $ U.slice (ix * maxD2) maxD2 tab
      where
        step !acc !i !num = acc - fromIntegral num * unPr p
          where
            (!d1,!d2) = i `quotRem` maxD
            p = g1 + Pr (- fromIntegral d1) * g2 + Pr (- fromIntegral (d1+d2)) * g3


-- | These are the top three mappable regions for each chromsome, taken
-- from Heng Li's "filt35_50".  They total just shy of 2MB altogether;
-- we'll use them for parameter fitting.

good_regions :: [( Bytes, [(Int,Int)] )]
good_regions =
    [ ("1", [( 87790801, 87822909),(  2975248,  3003748),( 49175564, 49202682)])
    , ("2", [(172933866,172969274),( 72016632, 72047739),(179421187,179450188)])
    , ("3", [(126700279,126741143),( 47025806, 47056328),( 49678222, 49708583)])
    , ("4", [( 20244708, 20281684),(111534379,111568539),(105402662,105430227)])
    , ("5", [(134654381,134687221),( 94195249, 94225025),( 66378246, 66405287)])
    , ("6", [( 50781199, 50816527),( 47822738, 47851832),( 93965080, 93992193)])
    , ("7", [( 50443055, 50479884),( 27127108, 27154444),( 50707378, 50734468)])
    , ("8", [(140633542,140663973),(143601826,143629417),(142668520,142695917)])
    , ("9", [(137210024,137260825),( 23816847, 23850736),(124513463,124539970)])
    , ("10",[(130630869,130657212),( 11191315, 11216898),(131751339,131776023)])
    , ("11",[(  2701423,  2730710),( 15625359, 15653856),( 17606351, 17634725)])
    , ("12",[(  5601165,  5634882),( 54342242, 54373433),( 65886788, 65913420)])
    , ("13",[( 44578644, 44605366),(112699565,112722141),( 67789099, 67811513)])
    , ("14",[( 22842926, 22872590),(101962136,101990258),( 33968191, 33992602)])
    , ("15",[( 36124467, 36148530),( 97239853, 97262974),( 70369320, 70392391)])
    , ("16",[( 65235334, 65263721),( 86516836, 86544522),(  1811705,  1830767)])
    , ("17",[( 80003858, 80024543),( 32778941, 32798278),(  7306011,  7325155)])
    , ("18",[( 72958402, 73005526),( 31218235, 31248142),( 76730933, 76757908)])
    , ("19",[( 42650207, 42668254),( 42781659, 42799671),( 31832574, 31850483)])
    , ("20",[( 23007246, 23035150),( 57188413, 57210981),( 61147962, 61169082)])
    , ("21",[( 16420477, 16442980),( 17158618, 17179367),( 39746588, 39766281)])
    , ("22",[( 23433574, 23463899),( 19954339, 19981883),( 19729058, 19753981)]) ]
    -- ("X", [(152736264,152762436),(147879564,147904235),( 68085990, 68110400)])
    -- ("Y", [( 22551581, 22561871),( 23158377, 23167802),( 22137261, 22146193)]) ]

-- | Parameter estimation for a single sample.  The parameters are
-- divergence and heterozygosity.  We tabulate the data here and do the
-- estimation afterwards.  Returns the product of the
-- parameter-independent parts of the likehoods and the histogram
-- indexed by D and H (see @genotyping.pdf@ for details).
tabulateSingle :: (Functor m, MonadIO m) => Iteratee [Calls] m DivTable
tabulateSingle = do
    tab <- liftIO $ M.replicate (12 * maxD * maxD) (0 :: Int)
    DivTable <$> foldStreamM (\acc -> accum tab acc . p_snp_pile) (0 :: Double)
             <*> liftIO (U.unsafeFreeze tab)
  where
    maxD = 64

    -- We need GL values for the invariant, the three homozygous variant
    -- and the three single-event heterozygous variant cases.  The
    -- ordering is like in BCF, with the reference first.
    -- Ref ~ A ==> PL ~ AA, AC, CC, AG, CG, GG, AT, CT, GT, TT
    {-# INLINE accum #-}
    accum !tab !acc (Snp_GLs !gls !ref)
        | U.length gls /= 10           = error "Ten GL values expected for SNP!"      -- should not happen
        | ref `elem` [nucsC,nucsG]     = accum' 0 tab acc gls ref
        | ref `elem` [nucsA,nucsT]     = accum' 6 tab acc gls ref
        | otherwise                    = return acc                                   -- unknown reference

    -- The simple 2D table didn't work, it lacked resolution in some
    -- cases.  We make six separate tables instead so we can store two
    -- differences with good resolution in every case.
    {-# INLINE accum' #-}
    accum' !refix !tab !acc !gls !ref
        | g_RR >= g_RA && g_RA >= g_AA = store 0 g_RR g_RA g_AA
        | g_RR >= g_AA && g_AA >= g_RA = store 1 g_RR g_AA g_RA
        | g_RA >= g_RR && g_RR >= g_AA = store 2 g_RA g_RR g_AA
        | g_RA >= g_AA && g_AA >= g_RR = store 3 g_RA g_AA g_RR
        | g_RR >= g_RA                 = store 4 g_AA g_RR g_RA
        | otherwise                    = store 5 g_AA g_RA g_RR

      where
        g_RR | ref == nucsT = unPr $  U.unsafeIndex gls 9
             | ref == nucsG = unPr $  U.unsafeIndex gls 5
             | ref == nucsC = unPr $  U.unsafeIndex gls 2
             | otherwise    = unPr $  U.unsafeIndex gls 0

        g_RA                = unPr $ (U.unsafeIndex gls 1 + U.unsafeIndex gls 3 + U.unsafeIndex gls 6) / 3

        g_AA | ref == nucsT = unPr $ (U.unsafeIndex gls 0 + U.unsafeIndex gls 2 + U.unsafeIndex gls 5) / 3
             | ref == nucsG = unPr $ (U.unsafeIndex gls 0 + U.unsafeIndex gls 2 + U.unsafeIndex gls 9) / 3
             | ref == nucsC = unPr $ (U.unsafeIndex gls 0 + U.unsafeIndex gls 5 + U.unsafeIndex gls 9) / 3
             | otherwise    = unPr $ (U.unsafeIndex gls 2 + U.unsafeIndex gls 5 + U.unsafeIndex gls 9) / 3

        store !t !a !b !c = do let d1 = min (maxD-1) . round $ a - b
                                   d2 = min (maxD-1) . round $ b - c
                                   ix = (t + refix) * maxD * maxD + d1 * maxD + d2
                               liftIO $ M.read tab ix >>= M.write tab ix . succ
                               return $! acc + a

