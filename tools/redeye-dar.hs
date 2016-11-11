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

import Bio.Adna
import Bio.Bam
import Bio.Bam.Pileup
import Bio.Genocall
import Bio.Genocall.Estimators
import Bio.Prelude
import Bio.Util.AD
import Bio.Util.Pretty
import System.Console.GetOpt

import qualified Data.HashMap.Strict as H
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed.Mutable as U

data Conf = Conf {
    conf_report :: String -> IO (),
    conf_params :: Parameters }

defaultConf :: Conf
defaultConf = Conf 30 (\_ -> return ()) quietParameters

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "v"  ["verbose"]      (NoArg      set_verbose) "Print progress reports",
    Option "h?" ["help","usage"] (NoArg       disp_usage) "Print this message and exit" ]
  where
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

    model0 <- newIORef H.empty
    -- Meh.  New concept.  We'll operate repeatedly on a small set of
    -- regions.  For the time being, that region will have to be put
    -- into a suitable file, so here we expect precisely one file.  Read
    -- group annotations have to be present, we will apply one damage
    -- model per read group.  Have to decide when we're done.

    -- For each iteration:  read the input, decompose, pileup.  The
    -- "prior" damage model is the usual 'SubstModel', the "posterior"
    -- damage model needs to be a mutable 'MSubstModel'.  We feed
    -- likelihoods into the div/het estimation (just as in
    -- 'redeye-pileup'), but also into a caller that will estimate
    -- damage.
    tab <- emIter model0 conf_report files
    estimateSingle tab >>= \(u,v) -> pprint u >> pprint v


emIter :: IORef (HashMap Bytes (SubstModel, MSubstModel))
       -> (String -> IO ()) -> [FilePath] -> IO DivTable
emIter model report infiles =
        fmap fst                                                                    $
        concatInputs infiles >=> run                                                $ \hdr ->
        concatMapStreamM (decompose_dmg_from model)                                =$
        progressPos (\(rs, p, _) -> (rs, p)) "GT call at " report (meta_refs hdr)  =$
        pileup                                                                     =$
        mapStream calls                                                            =$
        zipStreams tabulateSingle skipToEof


calls :: Pile ( Mat44D, b ) -> Calls
calls pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !s = simple_snp_call   2 $ map (second (fmap fst)) $ uncurry (++) $ p_snp_pile pile
    !i = simple_indel_call 2 $ map (second (second (map (fmap fst)))) $ p_indel_pile pile


{-# INLINE decompose_dmg_from #-}
decompose_dmg_from :: IORef (HashMap Bytes (SubstModel, MSubstModel)) -> BamRaw -> IO [ PosPrimChunks (Mat44D, MMat44D) ]
decompose_dmg_from ref raw = do
    hm <- readIORef ref
    let rg = extAsString "RG" (unpackBam raw)
    model <- case H.lookup rg hm of
                Just mm -> return mm
                Nothing -> do mm <- (,) model0 <$> freeze_subst_model model0
                              writeIORef ref $! H.insert rg mm hm
                              return mm
    return $ maybe [] (:[]) $ decompose (from_mmodel model) raw

  where
    from_mmodel (m,mm) i
        | i >= 0 &&   i  <  V.length  (left_substs m) = ( V.unsafeIndex (left_substs   m)  i
                                                        , V.unsafeIndex (left_substs  mm)  i )
        | i <  0 && (-i) >= V.length (right_substs m) = ( V.unsafeIndex (right_substs  m) (-i-1)
                                                        , V.unsafeIndex (right_substs mm) (-i-1) )
        | otherwise                                   = ( middle_substs m, middle_substs mm )

    idmat   = scalarMat 1

    model0 = SubstModel { left_substs   = V.replicate 12 idmat
                        , middle_substs =                idmat
                        , right_substs  = V.replicate 12 idmat }

    nullmat = MMat44D <$> U.replicate 16 (0::Double)
    freeze_subst_model m = SubstModel <$> V.mapM (const nullmat) (left_substs   m)
                                      <*>               nullmat
                                      <*> V.mapM (const nullmat) (right_substs  m)

