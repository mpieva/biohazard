{-# LANGUAGE RecordWildCards, NamedFieldPuns #-}
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

import Bio.Bam
import Bio.Genocall.Damage
import Bio.Prelude
import Bio.Util.Pretty
import Bio.Util.AD
import System.Console.GetOpt

data Conf = Conf {
    conf_lmin :: Int,
    conf_report :: String -> IO (),
    conf_params :: Parameters }

defaultConf :: Conf
defaultConf = Conf 25 (\_ -> return ()) quietParameters

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "m"  ["min-length"]   (ReqArg set_lmin  "LEN") "Set minimum length to LEN (25)",
    Option "v"  ["verbose"]      (NoArg      set_verbose) "Print progress reports",
    Option "h?" ["help","usage"] (NoArg       disp_usage) "Print this message and exit" ]
  where
    set_lmin   a c = readIO a >>= \l -> return $ c { conf_lmin = l }
    set_verbose  c = return $ c { conf_report = hPutStrLn stderr, conf_params = debugParameters }

    disp_usage  _ = do pn <- getProgName
                       let blah = "Usage: " ++ pn ++ " [OPTION...] [LIBRARY-NAME...]"
                       putStrLn $ usageInfo blah options
                       exitSuccess

main :: IO ()
main = do
    (opts, files, errors) <- getOpt Permute options <$> getArgs
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Conf{..} <- foldl (>>=) (return defaultConf) opts

    -- XXX  We currently don't get confidence intervals.  Are they even
    -- all that interesting?
    estimateDamageFromFiles conf_lmin conf_params files >>= pprint

    -- XXX  Many optimizations fit only one parameter.  Newton-Iteration
    -- should be more efficient that the generic CG method.
