{-# LANGUAGE RecordWildCards, BangPatterns, OverloadedStrings, FlexibleContexts #-}
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
-- overhead of 150% useless likelihood values for the sex chromosomes.

-- XXX Narf... for the parameter estimate, we need the reference.
-- Either supply it, or derive it from the MD field.  Except pileup
-- doesn't implement that.  :(

import Bio.Base
import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Bam.Pileup
import Bio.Genocall
import Bio.Genocall.Adna
import Bio.Genocall.AvroFile
import Bio.Iteratee
import Control.Applicative
import Control.Monad
import Data.Avro
import System.Console.GetOpt
import System.Environment
import System.Exit
import System.IO

import qualified Data.ByteString.Char8          as S
import qualified Data.Iteratee                  as I
import qualified Data.Text.Encoding             as T
import qualified Data.Vector.Unboxed            as U

type OIter = Conf -> Refs -> Iteratee [Calls] IO ()

data Conf = Conf {
    conf_output      :: (OIter -> IO ()) -> IO (),
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

    conf_output $ \oiter ->
        mergeInputs combineCoordinates files >=> run $ \hdr ->
            filterStream (not . br_isUnmapped)                      =$
            filterStream (isValidRefseq . br_rname)                 =$
            progressPos "GT call at " conf_report (meta_refs hdr)      =$
            pileup dmg_model                                        =$
            by_groups p_refseq (\rs out -> do
                let !sname = sq_name $ getRef (meta_refs hdr) rs
                liftIO $ conf_report $ S.unpack sname
                mapStream (calls conf_theta) out)                   =$
            oiter conf (meta_refs hdr)


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
    !s = (U.zipWith (*) x y, if r == nucsN then r' else r)
      where
        (x,r)  = (maq_snp_call 2 theta $ fst $ p_snp_pile pile)
        (y,r') = (maq_snp_call 2 theta $ snd $ p_snp_pile pile)


by_groups :: ( Monad m, ListLike s a, Nullable s, Eq b ) => (a -> b) -> (b -> Enumeratee s t m r) -> Enumeratee s t m r
by_groups f k out = do
    mhd <- peekStream
    case fmap f $ mhd of
        Nothing -> return out
        Just b0 -> takeWhileE ((==) b0 . f) =$ k b0 out >>= by_groups f k

output_avro :: Handle -> Conf -> Refs -> Iteratee [Calls] IO ()
output_avro hdl _cfg refs = compileBlocks refs =$
                            writeAvroContainer ContainerOpts{..} =$
                            mapChunksM_ (S.hPut hdl)
  where
    objects_per_block = 16      -- XXX should be more?
    filetype_label = "Genotype Likelihoods V0.1"


-- Serialize the results from genotype calling in a sensible way.  We
-- write an Avro file, but we add another blocking layer on top so we
-- don't need to endlessly repeat coordinates.

compileBlocks :: Monad m => Refs -> Enumeratee [Calls] [GenoCallBlock] m a
compileBlocks refs = convStream $ do
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
                    { reference_name = T.decodeLatin1 $ sq_name $ getRef refs rs
                    , start_position = p0
                    , called_sites   = reverse acc } ]

    pack c1 = rlist indel_variants `seq` GenoCallSite{..}
      where
        !snp_stats         = p_snp_stat c1
        !indel_stats       = p_indel_stat c1
        !snp_likelihoods   = compact_likelihoods $ fst $ p_snp_pile c1
        !indel_likelihoods = compact_likelihoods $ fst $ p_indel_pile c1
        !ref_allele        = snd $ p_snp_pile c1
        !indel_variants    = snd $ p_indel_pile c1

        rlist [] = ()
        rlist (x:xs) = x `seq` rlist xs

