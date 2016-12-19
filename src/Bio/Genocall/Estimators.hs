{-# LANGUAGE DeriveGeneric #-}
{-# OPTIONS_GHC -O0 #-}
module Bio.Genocall.Estimators (
        tabulateSingle,
        estimateSingle,
        DivEst(..),
        ExtModel(..),
        DivTable(..),
        good_regions
    ) where

-- XXX  Lots of (currently) dead code that needs to be removed.

-- XXX  Many optimizations fit only one parameter.  Newton-Iteration
-- should be more efficient than the generic CG method.
--
-- For Newton-Iteration, the update is 'x := x + dt' where
-- 'dt = -H(x)^{-1} \Nabla f(x)' or 'H(x) dt = -\Nabla f(x)'

-- import Bio.Adna
import Bio.Bam
import Bio.Bam.Pileup                ( p_snp_pile )
import Bio.Genocall                --  ( Snp_GLs(..), Calls, SubstModel )
import Bio.Prelude
import Bio.Util.AD
import Bio.Util.AD2
import Bio.Util.Numeric              ( log1pexp )
import Data.Aeson
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

{- smat :: (Floating a, IsDouble a) => Double -> a -> a -> a -> a -> [a]
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
    pT = fromDouble $ 0.5 - gc * 0.5 -}


-- | Estimates GC fraction.  Base composition is stable with respect to
-- our substitution model, so we should be able to estimate it from the
-- reference sequence.  Practically, since G and C have to be balanced,
-- we just estimate the GC fraction.

{- est_gc_frac :: DmgStats a -> Double
est_gc_frac DmgStats{..} = fromIntegral gc_counts / fromIntegral (at_counts+gc_counts)
  where
    getV n = maybe 0 U.sum . lookup (Just n)
    at_counts = sum [ getV n v | n <- [nucA,nucT], v <- [basecompo5, basecompo3] ]
    gc_counts = sum [ getV n v | n <- [nucG,nucC], v <- [basecompo5, basecompo3] ] -}

-- | Estimates mutation rates for transitions and transversion.  We
-- write down the likelihood function with just two general deamination
-- rates.  Likewise, our substitution statistics collapse to one single
-- matrix.  This is not perfect, but hopefully good enough to get a
-- quick estimate of the two mutation rates.  (Numeric optimization,
-- sorry it's in IO.)
{- est_mut_rate :: Parameters -> Double -> DmgStats a -> IO (Double, Double)
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
    lk_fun _ = error "Inconceivable!"  -}


-- | Estimate deamination rates.  We take the base composition (GC
-- fraction) and the two mutation rates (ti,tv) as constant, then
-- estimate two deamination rates (C->T and G->A) for a given number of
-- positions, and then the remainder.  (They are now decoupled, so we
-- could as well estimate the parameters (just two!) using Newton
-- iteration.)
{- est_deamination :: Parameters -> Int -> Double -> (Double,Double) -> DmgStats a
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
    lk_fun _ _ = error "Inconceivable!" -}

data DivTable = DivTable !Double !(U.Vector Int) deriving (Show, Generic)

instance Binary DivTable

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

instance ToJSON   DivEst
instance FromJSON DivEst

data ExtModel = ExtModel { population :: DivEst
                         , pop_separate :: Maybe DivEst
                         , damage :: SubstModels }
  deriving (Show, Generic)

instance ToJSON ExtModel
instance FromJSON ExtModel

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


-- | Parameter estimation for a single sample.  The parameters are
-- divergence and heterozygosity.  We tabulate the data here and do the
-- estimation afterwards.  Returns the product of the
-- parameter-independent parts of the likehoods and the histogram
-- indexed by D and H (see @genotyping.pdf@ for details).
tabulateSingle :: MonadIO m => Iteratee [Calls] m DivTable
tabulateSingle = do
    tab <- liftIO $ M.replicate (12 * maxD * maxD) (0 :: Int)
    DivTable `liftM` foldStreamM (\acc -> accum tab acc . p_snp_pile) (0 :: Double)
                `ap` liftIO (U.unsafeFreeze tab)
  where
    maxD = 64

    -- We need GL values for the invariant, the three homozygous variant
    -- and the three single-event heterozygous variant cases.  The
    -- ordering is like in BCF, with the reference first.
    -- Ref ~ A ==> PL ~ AA, AC, CC, AG, CG, GG, AT, CT, GT, TT
    {-# INLINE accum #-}
    accum !tab !acc !snp
        | ref `elem` [nucsC,nucsG]     = accum' 0 tab acc gls ref
        | ref `elem` [nucsA,nucsT]     = accum' 6 tab acc gls ref
        | otherwise                    = return acc                                   -- unknown reference
      where
        ref = snp_refbase snp
        gls = snp_gls     snp

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

-- | These are the longest contiguous mappable regions (chromsome, start
-- position, length) on the autosomes taken from Heng Li's "filt35_50"
-- until their cumulative length is greater than 10MB.  We'll use them
-- for parameter fitting.

good_regions :: [( Bytes, Int, Int )]
good_regions = [ ("9",137210024,50801), ("18",72958402,47124), ("3",126700279,40864), ("4",20244708,36976),
    ("7",50443055,36829), ("2",172933866,35408), ("6",50781199,35328), ("4",111534379,34160), ("9",23816847,33889),
    ("12",5601165,33717), ("5",134654381,32840), ("1",87790801,32108), ("12",54342242,31191), ("2",72016632,31107),
    ("3",47025806,30522), ("8",140633542,30431), ("3",49678222,30361), ("22",23433574,30325), ("18",31218235,29907),
    ("5",94195249,29776), ("14",22842926,29664), ("3",10955664,29450), ("11",2701423,29287), ("6",47822738,29094),
    ("2",179421187,29001), ("1",2975248,28500), ("11",15625359,28497), ("16",65235334,28387), ("11",17606351,28374),
    ("11",6627401,28349), ("14",101962136,28122), ("11",131331738,28031), ("3",35676569,27919), ("20",23007246,27904),
    ("16",86516836,27686), ("8",143601826,27591), ("4",105402662,27565), ("22",19954339,27544), ("8",142668520,27397),
    ("7",27127108,27336), ("3",135073061,27276), ("1",49175564,27118), ("6",93965080,27113), ("7",50707378,27090),
    ("5",66378246,27041), ("18",76730933,26975), ("13",44578644,26722), ("12",65886788,26632), ("9",124513463,26507),
    ("11",31654911,26355), ("10",130630869,26343), ("3",96522468,26316), ("5",71281187,26300), ("18",72601061,26219),
    ("8",89116777,26217), ("8",77595904,26119), ("7",36636342,26036), ("2",72353320,25820), ("3",112988383,25762),
    ("7",155239629,25756), ("1",2789639,25749), ("1",48213405,25681), ("4",136124313,25649), ("1",196271156,25641),
    ("10",11191315,25583), ("2",206608755,25468), ("12",54410147,25352), ("2",12854107,25292), ("11",11581616,25288),
    ("12",13700919,25192), ("8",142864620,25147), ("2",179465162,25139), ("9",14514334,25107), ("3",35717320,25063),
    ("22",19729058,24923), ("9",109680443,24899), ("6",94115565,24860), ("8",93505503,24846), ("4",158790023,24704),
    ("10",131751339,24684), ("1",181704049,24484), ("10",80682212,24465), ("3",169164284,24463), ("11",17511342,24437),
    ("14",33968191,24411), ("1",53992234,24288), ("4",183056381,24222), ("2",71820551,24155), ("9",137265761,24075),
    ("15",36124467,24063), ("14",46406527,24053), ("1",83210498,24021), ("6",108476739,23964), ("5",131587752,23956),
    ("3",48615599,23886), ("6",55728531,23880), ("7",127697189,23869), ("12",50346181,23863), ("3",144228991,23831),
    ("4",8265418,23682), ("2",45148758,23679), ("5",134540523,23644), ("5",134371610,23626), ("1",42036301,23604),
    ("8",116656966,23590), ("9",96702236,23589), ("1",37307947,23570), ("1",210850735,23541), ("9",135459556,23527),
    ("2",51135518,23503), ("11",45670249,23503), ("3",122993389,23483), ("14",56745421,23461), ("12",5898019,23455),
    ("2",163037779,23446), ("11",69447933,23373), ("2",212035036,23357), ("6",96460435,23294), ("18",44764982,23243),
    ("10",101278523,23242), ("5",160036671,23241), ("1",49124538,23176), ("7",18260998,23158), ("2",99167682,23147),
    ("14",29232562,23127), ("15",97239853,23121), ("10",50494407,23089), ("15",70369320,23071), ("3",168827492,23039),
    ("8",77761638,23015), ("15",38157092,22957), ("10",11365418,22949), ("3",105255227,22893), ("1",34814158,22881),
    ("1",37392280,22867), ("11",71705177,22832), ("12",2418673,22823), ("1",72511401,22645), ("7",1265072,22602),
    ("13",112699565,22576), ("20",57188413,22568), ("21",16420477,22503), ("3",129287740,22468), ("9",125868895,22467),
    ("1",160050031,22456), ("12",48128061,22430), ("18",5282725,22418), ("13",67789099,22414), ("6",40287102,22408),
    ("18",5880352,22384), ("8",143528138,22367), ("12",85676068,22365), ("8",84320057,22352), ("9",37321783,22293),
    ("1",209775949,22272), ("3",18389755,22252), ("4",147944592,22243), ("8",133906982,22224), ("11",133329130,22217),
    ("3",173313782,22215), ("5",11889097,22179), ("11",79649083,22174), ("5",132731240,22137), ("7",157466401,22137),
    ("2",179634275,22116), ("11",43589949,22100), ("4",52766804,22099), ("12",3355258,22074), ("5",131397743,22070),
    ("10",80291692,21987), ("9",14296528,21964), ("10",98269731,21962), ("3",16914253,21948), ("1",197561249,21908),
    ("6",120906102,21884), ("2",193650467,21851), ("2",128384358,21847), ("3",48662757,21839), ("6",43733826,21830),
    ("4",10215804,21821), ("22",46353587,21773), ("1",239171473,21764), ("12",54373478,21741), ("15",88509335,21732),
    ("18",35145944,21709), ("15",95071676,21703), ("1",3220853,21688), ("1",82261961,21677), ("13",26246199,21672),
    ("2",121320633,21663), ("5",95664126,21658), ("18",53770513,21624), ("5",108573486,21615), ("10",1762508,21612),
    ("8",131982804,21550), ("12",59791869,21536), ("2",164194661,21524), ("11",64389512,21519), ("10",48420962,21516),
    ("12",55408280,21514), ("14",33881612,21508), ("8",93616886,21507), ("3",134838560,21418), ("1",32215568,21396),
    ("1",88245563,21394), ("5",163625223,21390), ("1",154971490,21367), ("1",2393882,21354), ("2",119583947,21343),
    ("10",125760979,21336), ("12",1928169,21298), ("11",15787032,21257), ("1",160020996,21256), ("14",34162064,21249),
    ("5",152857671,21234), ("10",132906774,21228), ("12",91441003,21228), ("1",168471368,21187), ("11",6660346,21183),
    ("15",93125306,21146), ("2",159569436,21136), ("7",114286123,21134), ("3",55497586,21131), ("20",61147962,21120),
    ("18",31314083,21087), ("5",87948855,21077), ("12",107266456,21071), ("1",112382584,21045), ("1",44871084,21041),
    ("8",59833303,21037), ("20",12185890,21016), ("1",175835258,20982), ("11",132197740,20962), ("6",128317896,20960),
    ("1",216882694,20949), ("11",61329536,20948), ("4",24007933,20895), ("7",42485962,20843), ("10",44446025,20828),
    ("11",2153220,20823), ("3",52112801,20818), ("1",75232624,20802), ("1",81018216,20782), ("21",17158618,20749),
    ("11",15830777,20738), ("20",11895285,20734), ("10",52998525,20713), ("5",93635203,20708), ("17",80003858,20685),
    ("9",129968922,20683), ("6",70561018,20671), ("4",97130238,20662), ("8",135136261,20658), ("2",23903685,20632),
    ("22",24548899,20628), ("9",7497211,20620), ("13",110428506,20601), ("15",48384378,20563), ("12",70712242,20555),
    ("7",149504838,20520), ("1",42189196,20518), ("2",109987828,20491), ("3",18073619,20438), ("1",239875096,20430),
    ("11",124937110,20416), ("5",11382753,20401), ("20",39447354,20400), ("14",105252979,20394), ("9",24026887,20386),
    ("5",132623090,20381), ("8",65608184,20365), ("6",33743378,20351), ("2",151329557,20344), ("20",37892360,20340),
    ("1",3367177,20331), ("5",122512407,20331), ("10",102974441,20331), ("13",112755319,20326), ("3",12992872,20315),
    ("13",42170961,20306), ("10",80874296,20300), ("20",60444066,20292), ("2",73146343,20284), ("3",50190348,20284),
    ("22",50340074,20277), ("1",77030540,20253), ("13",102679591,20222), ("10",81040310,20197), ("9",9211551,20188),
    ("18",22738887,20181), ("3",44052705,20176), ("2",54880032,20166), ("5",37822780,20150), ("18",53007829,20135),
    ("9",129187481,20128), ("6",144270559,20091), ("1",205403570,20090), ("4",44421876,20045), ("10",44861946,20023),
    ("4",67049100,19998), ("7",32096011,19967), ("7",34131372,19967), ("1",110017253,19938), ("1",77080436,19922),
    ("5",88010450,19920), ("7",15716479,19911), ("15",33826577,19904), ("9",126521960,19877), ("14",58022067,19853),
    ("7",27160875,19845), ("4",90823485,19825), ("18",45240281,19818), ("1",3005635,19809), ("2",157007834,19789),
    ("12",114832107,19789), ("6",45288303,19788), ("6",99074944,19783), ("11",2670228,19779), ("5",73567689,19774),
    ("6",51224143,19774), ("5",102223202,19766), ("15",89899715,19750), ("4",104970277,19740), ("15",70349377,19722),
    ("6",101837270,19718), ("4",96325943,19709), ("2",168098209,19708), ("21",39746588,19693), ("3",59949020,19686),
    ("2",119358552,19681), ("2",67497291,19679), ("12",105117749,19661), ("3",180017017,19638), ("4",1793613,19638),
    ("7",92453277,19635), ("20",16315532,19598), ("5",169493896,19574), ("3",18469949,19564), ("4",2061136,19548),
    ("1",175528763,19535), ("15",87190931,19514), ("14",101991199,19498), ("9",9497024,19493), ("11",15896925,19493),
    ("22",43629376,19444), ("11",64363871,19436), ("3",50399986,19433), ("1",244210029,19431), ("11",2591551,19424),
    ("11",1239413,19417), ("9",116776003,19416), ("1",198338602,19414), ("7",8501937,19411), ("1",120452998,19403),
    ("14",24722981,19395), ("15",39870070,19395), ("2",134164180,19386), ("11",109293469,19365), ("3",113929897,19362),
    ("4",107835785,19338), ("17",32778941,19337), ("5",122424884,19327), ("13",105234793,19317), ("21",18170750,19309),
    ("5",134719430,19301), ("14",56588948,19294), ("6",39449811,19289), ("7",116762432,19282), ("8",79510243,19281),
    ("11",75301993,19266), ("3",89452118,19229), ("6",40123796,19208), ("1",181452626,19207), ("14",33860687,19197),
    ("11",124731491,19188), ("8",62768789,19187), ("9",124458437,19180), ("9",4097024,19172), ("2",205522223,19160),
    ("1",166027772,19156), ("10",130504847,19155), ("17",7306011,19144), ("10",23479527,19123), ("4",2245744,19117),
    ("3",46901850,19109), ("3",139771287,19105), ("9",34548966,19104), ("8",110811609,19093), ("2",220339781,19077),
    ("16",1811705,19062), ("7",155530076,19061), ("3",110181036,19060), ("2",127805324,19050), ("7",55182828,19047),
    ("7",27223390,19039), ("1",51432006,19038), ("4",88753129,19035), ("15",63122126,19030), ("5",87837753,19018),
    ("2",36665285,19015), ("18",35186913,18942), ("5",129172415,18940), ("15",68111506,18935), ("7",25440487,18925),
    ("8",93889992,18924), ("17",56398946,18921), ("13",107170684,18912), ("13",102563238,18908), ("10",43743824,18896),
    ("15",50211905,18888), ("10",72530467,18883), ("11",113271284,18878), ("12",81086307,18874), ("1",175352300,18867),
    ("4",10071903,18866), ("14",89644782,18857), ("3",50419448,18855), ("11",93212299,18848), ("7",44268521,18838),
    ("12",123458577,18812), ("5",61089957,18811), ("2",19554179,18803), ("7",132105968,18781), ("15",97129465,18776),
    ("20",21492872,18766), ("2",149670253,18746), ("10",49649785,18737), ("10",125062922,18725), ("22",51154624,18693),
    ("20",57815879,18683), ("2",76497569,18675), ("8",96701949,18673), ("1",87612820,18670), ("5",158974585,18664),
    ("13",58980956,18664), ("2",156229090,18654), ("5",159334986,18638), ("14",101316928,18631), ("3",42692075,18627),
    ("7",55163429,18609), ("17",63881123,18590), ("9",73721557,18579), ("1",177235920,18568), ("3",10800230,18565),
    ("2",66495731,18555), ("2",6162803,18553), ("6",102243030,18543), ("5",160945102,18539), ("7",94527596,18539),
    ("11",128360895,18538), ("5",113693991,18537), ("3",127979662,18534) ]
