{-# LANGUAGE BangPatterns, RecordWildCards, OverloadedStrings #-}
-- Here we read the tables from sample_div_tables, add them up as
-- necessary, estimate divergence and heterozygosity from them, and
-- store the result back.  The estimate can be done for regions, which
-- are defined by regular expressions.

import Bio.Prelude
import Bio.Genocall.Estimators
import Bio.Util.AD
import Bio.Util.AD2
import Bio.Util.Numeric              ( log1p )
import Bio.Util.Regex                ( regComp, regMatch )
import Control.Concurrent.Async      ( async, wait )
import System.Console.GetOpt

import qualified Data.HashMap.Strict            as H
import qualified Data.Vector.Storable           as VS
import qualified Data.Vector.Unboxed            as U

data Conf = Conf { conf_metadata :: FilePath
                 , conf_regions  :: [Text]
                 , conf_verbose  :: Bool }

defaultConf :: Conf
defaultConf = Conf (error "no metadata file specified") [] False

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "c"  ["config"] (ReqArg set_conf    "FILE") "Set name of json config file to FILE",
    Option "r"  ["region"] (ReqArg add_region "REGEX") "What matches REGEX becomes a region",
    Option "H"  ["human"]            (NoArg set_human) "Use regions for a human genome",
    Option "v"  ["verbose"]         (NoArg be_verbose) "Print more diagnostics",
    Option "h?" ["help","usage"]    (NoArg disp_usage) "Display this message" ]
  where
    be_verbose       c = return $ c { conf_verbose = True }
    set_conf      fn c = return $ c { conf_metadata = fn }
    add_region    re c = return $ c { conf_regions = fromString re : conf_regions c }
    set_human        c = return $ c { conf_regions = [ "^(chr)?[0-9]+$", "^(chr)?X$", "^(chr)?Y$" ] }

    disp_usage _ = do pn <- getProgName
                      let blah = "Usage: " ++ pn ++ " [OPTION...] [SAMPLE...]"
                      putStrLn $ usageInfo blah options
                      exitSuccess


main :: IO ()
main = do
    (opts, samples, errs) <- getOpt Permute options `fmap` getArgs
    Conf{..} <- foldl (>>=) (return defaultConf) opts
    unless (null errs) $ mapM_ (hPutStrLn stderr) errs >> exitFailure

    meta0 <- readMetadata conf_metadata

    let eff_samples = if null      samples then H.keys meta0     else map fromString samples
        eff_regions = if null conf_regions then [""]             else conf_regions

    updates <- forM eff_samples >=> mapM wait . concat $ \sample -> case H.lookup sample meta0 of

            Nothing -> do hPutStrLn stderr $ "unknown sample " ++ show sample
                          exitFailure

            Just smp -> forM eff_regions $ \rgn -> async $ do
                            (!v,!m) <- estimateSingle
                                       $ foldl1' (\(DivTable a u) (DivTable b v) -> DivTable (a+b) (U.zipWith (+) u v))
                                       $ H.elems
                                       $ H.filterWithKey (match rgn)
                                       $ sample_div_tables smp

                            when conf_verbose $ putStrLn $
                                "Estimate done for " ++ show sample ++ ", " ++ show rgn ++ ":\n" ++ unpack m
                            return (rgn, sample, v)

    let upd1 mdata (sample, rgn, val) = H.adjust upd_smp sample mdata
          where upd_smp smp' = smp' { sample_divergences = H.insert rgn val $ sample_divergences smp' }

    updateMetadata (\mdata -> foldl' upd1 mdata updates) conf_metadata

  where
    match :: Text -> Text -> a -> Bool
    match rgn = const . regMatch (regComp $ unpack rgn) . unpack

