{-# LANGUAGE RecordWildCards, BangPatterns, OverloadedStrings, FlexibleContexts, RankNTypes #-}
-- Command line driver for simple genotype calling.  We have three
-- separate steps:  Pileup from a BAM file (or multiple merged files) to
-- produce likelihoods (and some auxillary statistics).  These are
-- written into an Avro container.  Next we need to estimate parameters,
-- in the simplest case divergence and heterozygosity.  We can save some
-- time by fusing this with the first step.  The final step is calling
-- bases by scnaning the Avro container and applying some model, and
-- again, in the simplest case that's just divergence and
-- heterozygosity.  We keep that separate, because different models will
-- require different programs.  So here we produce likelihoods and
-- a simple model fit.

-- The likelihoods depend on damage parameters and an error model,
-- otherwise they are 'eternal'.  (For the time being, it's probably
-- wise to go with the naïve error model.)  Technically, they also
-- depend on ploidy, but since only diploid organisms are interesting
-- right now, we fix that to two.  We pay some overhead on the sex
-- chromosomes, but the simplification is worth it.

-- About damage parameters:  We effectively have three different models
-- (SS, DS, no damage) and it may not be possible to choose one a
-- priori.  To manage this cleanly, we should have one universal model,
-- but the three we have are not generalizations of each other.
-- However, all can be generalized into one model with slightly more
-- parameters.  See tools/dmg-est.hs for how we fit the model.

-- Calling is always diploid, for maximum flexibility.  We don't really
-- support higher ploidies, so the worst damage is that we output an
-- overhead of 150% useless likelihood values for the sex chromosomes
-- and maybe estimate heterozygosity where there is none.

import AD
import Bio.Base
import Bio.Bam.Header
import Bio.Bam.Pileup
import Bio.Bam.Reader
import Bio.Bam.Rec
import Bio.Genocall
import Bio.Genocall.Adna
import Bio.Genocall.AvroFile
import Bio.Iteratee
import Bio.Util                     ( log1p )
import Control.Applicative
import Control.Monad
import Data.Aeson
import Data.Avro
import Data.Ix
import System.Console.GetOpt
import System.Environment
import System.Exit
import System.IO

import qualified Data.ByteString.Char8          as S
import qualified Data.ByteString.Lazy           as BL
import qualified Data.Foldable                  as F
import qualified Data.HashMap.Strict            as H
import qualified Data.Text.Encoding             as T
import qualified Data.Vector.Storable           as VS
import qualified Data.Vector                    as V
import qualified Data.Vector.Unboxed            as U
import qualified Data.Vector.Unboxed.Mutable    as M

type OIter = Conf -> Refs -> Iteratee [Calls] IO ()

data Conf = Conf {
    conf_output      :: forall r . (OIter -> IO r) -> IO r,
    conf_damage      :: Maybe (DamageParameters Double),
    conf_loverhang   :: Maybe Double,
    conf_roverhang   :: Maybe Double,
    conf_ds_deam     :: Double,
    conf_ss_deam     :: Double,
    conf_theta       :: Maybe Double,
    conf_report      :: String -> IO () }

defaultConf :: Conf
defaultConf = Conf ($ output_avro stdout) Nothing Nothing Nothing 0.02 0.45 Nothing (\_ -> return ())

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o" ["output"]        (ReqArg set_avro_out  "FILE") "Write AVRO output to FILE",
    Option "D" ["damage"]        (ReqArg set_damage   "PARMS") "Set universal damage parameters",
    Option "l" lover_param       (ReqArg set_loverhang "PROB") "Parameter for 5' overhang length is PROB",
    Option "r" rover_param       (ReqArg set_roverhang "PROB") "Parameter for 3' overhang length is PROB, assume single-strand prep",
    Option "d" ds_param          (ReqArg set_ds_deam   "FRAC") "Deamination rate in double stranded section is FRAC",
    Option "s" ss_param          (ReqArg set_ss_deam   "FRAC") "Deamination rate in single stranded section is FRAC",
    Option "t" dep_param         (ReqArg set_theta     "FRAC") "Set dependency coefficient to FRAC (\"N\" to turn off)",
    Option "v" ["verbose"]       (NoArg            be_verbose) "Print more diagnostics",
    Option "h?" ["help","usage"] (NoArg            disp_usage) "Display this message" ]
  where
    lover_param = ["overhang-param","left-overhang-param"]
    rover_param = ["right-overhang-param"]
    ds_param    = ["deamination-rate","ds-deamination-rate","double-strand-deamination-rate"]
    ss_param    = ["ss-deamination-rate","single-strand-deamination-rate"]
    dep_param   = ["theta","dependency-coefficient"]

    disp_usage _ = do pn <- getProgName
                      let blah = "Usage: " ++ pn ++ " [OPTION...] [BAM-FILE...]"
                      putStrLn $ usageInfo blah options
                      exitFailure

    be_verbose       c = return $ c { conf_report = hPutStrLn stderr }
    set_avro_out "-" c = return $ c { conf_output = \k -> k $ output_avro stdout }
    set_avro_out  fn c = return $ c { conf_output = \k -> withFile fn WriteMode $ k . output_avro }

    set_theta    "N" c = return $ c { conf_theta  = Nothing }
    set_theta      a c = (\t -> c { conf_theta       = Just   t }) <$> readIO a

    set_loverhang  a c = (\l -> c { conf_loverhang   = Just   l }) <$> readIO a
    set_roverhang  a c = (\l -> c { conf_roverhang   = Just   l }) <$> readIO a
    set_ss_deam    a c = (\r -> c { conf_ss_deam     =        r }) <$> readIO a
    set_ds_deam    a c = (\r -> c { conf_ds_deam     =        r }) <$> readIO a
    set_damage     a c = (\u -> c { conf_damage      = Just   u }) <$> readIO a

main :: IO ()
main = do
    (opts, files, errs) <- getOpt Permute options <$> getArgs
    unless (null errs) $ mapM_ (hPutStrLn stderr) errs >> exitFailure
    conf@Conf{..} <- foldl (>>=) (return defaultConf) opts

    let no_damage   = conf_report "using no damage model" >> return noDamage
        ss_damage p = conf_report ("using single strand damage model with " ++ show p) >> return (univDamage p)
        ds_damage p = conf_report ("using double strand damage model with " ++ show p) >> return (univDamage p)
        u_damage  p = conf_report ("using universal damage parameters " ++ show p) >> return (univDamage p)

    dmg_model <- case (conf_damage, conf_loverhang, conf_roverhang) of
            (Just u,        _, _) -> u_damage u
            (_, Nothing, Nothing) -> no_damage
            (_, Just pl, Nothing) -> ds_damage $ DP 0 0 0 0 conf_ss_deam conf_ds_deam pl
            (_, Nothing, Just pr) -> ss_damage $ DP conf_ss_deam conf_ds_deam pr pr 0 0 0
            (_, Just pl, Just pr) -> ss_damage $ DP conf_ss_deam conf_ds_deam pl pr 0 0 0

    (tab,()) <- conf_output $ \oiter ->
        mergeInputs combineCoordinates files >=> run $ \hdr ->
            filterStream (not . isUnmapped . unpackBam)             =$
            filterStream (isValidRefseq . b_rname . unpackBam)      =$
            progressPos "GT call at " conf_report (meta_refs hdr)   =$
            pileup dmg_model                                        =$
            mapStream (calls conf_theta)                            =$
            zipStreams tabulateSingle (oiter conf $ meta_refs hdr)

    uncurry estimateSingle tab


-- | Ploidy is hardcoded as two here.  Can be changed if the need
-- arises.
--
-- XXX  For the time being, forward and reverse piles get concatenated.
-- For the naive call, this doesn't matter.  For the MAQ call, it feels
-- more correct to treat them separately and multiply (add?) the results.

calls :: Maybe Double -> Pile -> Calls
calls Nothing pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !s = simple_snp_call 2 $ uncurry (++) $ p_snp_pile pile
    !i = simple_indel_call 2 $ p_indel_pile pile

calls (Just theta) pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !i = simple_indel_call 2 $ p_indel_pile pile

    -- This lumps the two strands together
    -- !s = maq_snp_call 2 theta $ uncurry (++) $ p_snp_pile pile -- XXX

    -- This treats them separately
    !s | r == r'    = Snp_GLs (U.zipWith (*) x y) r     -- same ref base (normal case): multiply
       | r == nucsN = Snp_GLs y r'                      -- forward ref is N, use backward call
       | otherwise  = Snp_GLs x r                       -- else use forward call (even if this is incorrect,
      where                                             -- there is nothing else we can do here)
        Snp_GLs x r  = maq_snp_call 2 theta $ fst $ p_snp_pile pile
        Snp_GLs y r' = maq_snp_call 2 theta $ snd $ p_snp_pile pile


-- | Serialize the results from genotype calling in a sensible way.  We
-- write an Avro file, but we add another blocking layer on top so we
-- don't need to endlessly repeat coordinates.

compileBlocks :: Monad m => Enumeratee [Calls] [GenoCallBlock] m a
compileBlocks = convStream $ do
        c1 <- headStream
        tailBlock (p_refseq c1) (p_pos c1) (p_pos c1) . (:[]) $! pack c1
  where
    tailBlock !rs !p0 !po acc = do
        mc <- peekStream
        case mc of
            Just c1 | rs == p_refseq c1 && po+1 == p_pos c1 && po - p0 < 496 -> do
                    _ <- headStream
                    tailBlock rs p0 (po+1) . (:acc) $! pack c1

            _ -> return [ GenoCallBlock
                    { reference_name = rs
                    , start_position = p0
                    , called_sites   = reverse acc } ]

    pack c1 = rlist indel_variants `seq` GenoCallSite{..}
      where
        Snp_GLs snp_pls !ref_allele = p_snp_pile c1

        !snp_stats         = p_snp_stat c1
        !indel_stats       = p_indel_stat c1
        !snp_likelihoods   = compact_likelihoods snp_pls
        !indel_likelihoods = compact_likelihoods $ fst $ p_indel_pile c1
        !indel_variants    = snd $ p_indel_pile c1

        rlist [] = ()
        rlist (x:xs) = x `seq` rlist xs


output_avro :: Handle -> Conf -> Refs -> Iteratee [Calls] IO ()
output_avro hdl _cfg refs = compileBlocks =$
                            writeAvroContainer ContainerOpts{..} =$
                            mapChunksM_ (S.hPut hdl)
  where
    objects_per_block = 16      -- XXX should be more?
    filetype_label = "Genotype Likelihoods V0.1"
    initial_schemas = H.singleton "Refseq" $
        object [ "type" .= String "enum"
               , "name" .= String "Refseq"
               , "symbols" .= Array
                    (V.fromList . map (String . T.decodeUtf8 . sq_name) $ F.toList refs) ]
    meta_info = H.singleton "biohazard.refseq_length" $
                S.concat $ BL.toChunks $ encode $ Array $ V.fromList
                [ Number (fromIntegral (sq_length s)) | s <- F.toList refs ]


maxD :: Int
maxD = 63

rangeDs :: ((Int,Int,Int), (Int,Int,Int))
rangeDs = ((0,0,0),(11,maxD,maxD))

-- | Parameter estimation for a single sample.  The parameters are
-- divergence and heterozygosity.  We tabulate the data here and do the
-- estimation afterwards.  Returns the product of the
-- parameter-independent parts of the likehoods and the histogram
-- indexed by D and H (see @genotyping.pdf@ for details).
tabulateSingle :: (Functor m, MonadIO m) => Iteratee [Calls] m (Double, U.Vector Int)
tabulateSingle = do
    tab <- liftIO $ M.replicate (rangeSize rangeDs) (0 :: Int)
    (,) <$> foldStreamM (\acc -> accum tab acc . p_snp_pile) (0 :: Double)
        <*> liftIO (U.unsafeFreeze tab)
  where
    -- We need GL values for the invariant, the three homozygous variant
    -- and the three single-event heterozygous variant cases.  The
    -- ordering is like in BCF, with the reference first.
    -- Ref ~ A ==> PL ~ AA, AC, CC, AG, CG, GG, AT, CT, GT, TT
    {-# INLINE accum #-}
    accum !tab !acc (Snp_GLs !gls !ref)
        | U.length gls /= 10                   = error "Ten GL values expected for SNP!"      -- should not happen
        | ref `elem` [nucsC,nucsG]             = accum' 0 tab acc gls
        | ref `elem` [nucsA,nucsT]             = accum' 6 tab acc gls
        | otherwise                            = return acc                                   -- unknown reference

    -- The simple 2D table didn't work, it lacked resolution in some
    -- cases.  We make six separate tables instead so we can store two
    -- differences with good resolution in every case.
    {-# INLINE accum' #-}
    accum' refix !tab !acc !gls
        | g_RR >= g_RA && g_RA >= g_AA = store 0 g_RR g_RA g_AA
        | g_RR >= g_AA && g_AA >= g_RA = store 1 g_RR g_AA g_RA
        | g_RA >= g_RR && g_RR >= g_AA = store 2 g_RA g_RR g_AA
        | g_RA >= g_AA && g_AA >= g_RR = store 3 g_RA g_AA g_RR
        | g_AA >= g_RR && g_RR >= g_RA = store 4 g_AA g_RR g_RA
        | g_AA >= g_RA && g_RA >= g_RR = store 5 g_AA g_RA g_RR

      where
        g_RR = unPr $  U.unsafeIndex gls 0
        g_RA = unPr $ (U.unsafeIndex gls 1 + U.unsafeIndex gls 3 + U.unsafeIndex gls 6) / 3
        g_AA = unPr $ (U.unsafeIndex gls 2 + U.unsafeIndex gls 5 + U.unsafeIndex gls 9) / 3

        store t a b c = do let d1 = min maxD . round $ a - b
                               d2 = min maxD . round $ b - c
                               ix = index rangeDs (t + refix,d1,d2)
                           liftIO $ M.read tab ix >>= M.write tab ix . succ
                           return $! acc + a

estimateSingle :: Double -> U.Vector Int -> IO ()
estimateSingle llk_rr tab = do
    putStrLn $ "Estimating divergence parameters..."
    hFlush stdout
    (fit, res, stats) <- minimize quietParameters 0.0001 llk (U.fromList [0,0,0])
    let [delta, eta, eta2] = VS.toList fit
        dv = exp delta / (1 + exp delta)
        ht = exp eta / (1 + exp eta) -- * dv
        ht' = exp eta2 / (1 + exp eta2) -- * dv
        lk = finalValue stats - llk_rr
    putStrLn $ show res ++ ", " ++ show stats
    putStrLn $ "Divergence      = " ++ show dv
    putStrLn $ "Heterozygosity  = " ++ show ht
    putStrLn $ "Heterozygosity' = " ++ show ht'
    putStrLn $ "Log-Likelihood  = " ++ shows lk " ("
            ++ shows (finalValue stats) " - " ++ shows llk_rr ")"
  where
    llk :: [AD] -> AD
    llk [delta,eta,eta2] = llk' 0 delta eta + llk' 6 delta eta2

    llk' base delta eta = sum [ fromIntegral num * unPr (w_RR + w_AA * Pr (fromIntegral (d1+d2)) + w_RA * Pr (fromIntegral d1))
                              | i@(_,d1,d2) <- range ((base+0,0,0), (base+0,maxD,maxD))         -- g_RR g_RA g_AA
                              , let num = tab `U.unsafeIndex` index rangeDs i ]
                        + sum [ fromIntegral num * unPr (w_RR + w_AA * Pr (fromIntegral d1) + w_RA * Pr (fromIntegral (d1+d2)))
                              | i@(_,d1,d2) <- range ((base+1,0,0), (base+1,maxD,maxD))         -- g_RR g_AA g_RA
                              , let num = tab `U.unsafeIndex` index rangeDs i ]
                        + sum [ fromIntegral num * unPr (w_RR * Pr (fromIntegral d1) + w_AA * Pr (fromIntegral (d1+d2)) + w_RA)
                              | i@(_,d1,d2) <- range ((base+2,0,0), (base+2,maxD,maxD))         --  g_RA g_RR g_AA
                              , let num = tab `U.unsafeIndex` index rangeDs i ]
                        + sum [ fromIntegral num * unPr (w_RR * Pr (fromIntegral (d1+d2)) + w_AA * Pr (fromIntegral d1) + w_RA)
                              | i@(_,d1,d2) <- range ((base+3,0,0), (base+3,maxD,maxD))         --  g_RA g_AA g_RR
                              , let num = tab `U.unsafeIndex` index rangeDs i ]
                        + sum [ fromIntegral num * unPr (w_RR * Pr (fromIntegral d1) + w_AA + w_RA * Pr (fromIntegral (d1+d2)))
                              | i@(_,d1,d2) <- range ((base+4,0,0), (base+4,maxD,maxD))         --  g_AA g_RR g_RA
                              , let num = tab `U.unsafeIndex` index rangeDs i ]
                        + sum [ fromIntegral num * unPr (w_RR * Pr (fromIntegral (d1+d2)) + w_AA + w_RA * Pr (fromIntegral d1))
                              | i@(_,d1,d2) <- range ((base+5,0,0), (base+5,maxD,maxD))         --  g_AA g_RA g_RR
                              , let num = tab `U.unsafeIndex` index rangeDs i ]
      where
        !w_RR =        1 / Pr (log1p (exp delta))
        !w_AA = Pr delta / Pr (log1p (exp delta)) *      1 / Pr (log1p (exp eta))
        !w_RA = Pr delta / Pr (log1p (exp delta)) * Pr eta / Pr (log1p (exp eta))

