{-# LANGUAGE BangPatterns, RecordWildCards, OverloadedStrings #-}
-- Here we read the tables from sample_div_tables, add them up as
-- necessary, estimate divergence and heterozygosity from them, and
-- store the result back.  The estimate can be done for regions, which
-- are defined by regular expressions.

import Bio.Prelude
import Bio.Genocall.Metadata
import Bio.Util.AD
import Bio.Util.AD2
import Bio.Util.Numeric              ( log1p )
import Bio.Util.Regex                ( regComp, regMatch )
import Control.Concurrent.Async      ( async, wait )
import Numeric.LinearAlgebra.HMatrix ( eigSH', (><), toRows, scale )
import System.Console.GetOpt

import qualified Data.HashMap.Strict            as H
import qualified Data.Vector.Storable           as VS
import qualified Data.Vector.Unboxed            as U

data Conf = Conf { conf_metadata :: FilePath
                 , conf_regions  :: [Text]
                 , conf_purge    :: Bool
                 , conf_verbose  :: Bool }

defaultConf :: Conf
defaultConf = Conf (error "no metadata file specified") [] False False

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "c"  ["config"] (ReqArg set_conf    "FILE") "Set name of json config file to FILE",
    Option "r"  ["region"] (ReqArg add_region "REGEX") "What matches REGEX becomes a region",
    Option "p"  ["purge"]             (NoArg do_purge) "Purge tables after use",
    Option "H"  ["human"]            (NoArg set_human) "Use regions for a human genome",
    Option "v"  ["verbose"]         (NoArg be_verbose) "Print more diagnostics",
    Option "h?" ["help","usage"]    (NoArg disp_usage) "Display this message" ]
  where
    be_verbose       c = return $ c { conf_verbose = True }
    do_purge         c = return $ c { conf_purge = True }
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

    updates <- forM eff_samples >=> mapM wait $ \sample -> case H.lookup sample meta0 of

            Nothing -> do hPutStrLn stderr $ "unknown sample " ++ show sample
                          async $ return id

            Just smp -> async $ do
                ests <- forM eff_regions >=> mapM wait $ \rgn -> async
                                $ fmap ((,) rgn)
                                $ estimateSingle conf_verbose
                                $ foldl1' (\(DivTable a u) (DivTable b v) -> DivTable (a+b) (U.zipWith (+) u v))
                                $ H.elems
                                $ H.filterWithKey (match rgn)
                                $ sample_div_tables smp

                let app_purge = if conf_purge then appEndo (foldMap purge eff_regions) else id
                    upd_smp smp' = smp' { sample_divergences = ins_many ests $ sample_divergences smp'
                                        , sample_div_tables  = app_purge     $ sample_div_tables smp' }

                when conf_verbose $ putStrLn $ "Estimate done for " ++ show sample ++ "."
                return $ H.adjust upd_smp sample
    updateMetadata (foldr (.) id updates) conf_metadata

  where
    match :: Text -> Text -> a -> Bool
    match rgn = const . regMatch (regComp $ unpack rgn) . unpack

    purge :: Text -> Endo (H.HashMap Text a)
    purge rgn = Endo $ H.filterWithKey ((.) not . match rgn)

    ins_many :: [(Text,v)] -> H.HashMap Text v -> H.HashMap Text v
    ins_many = flip $ foldr (uncurry H.insert)


-- XXX we should estimate an indel rate, to be appended as the fourth
-- result (but that needs different tables)
estimateSingle :: Bool -> DivTable -> IO DivEst
estimateSingle verbose (DivTable llk_rr tab) = do
    (fit, res, stats) <- minimize quietParameters 0.0001 (llk tab) (U.fromList [0,0,0])
    let xform = map (\x -> exp x / (1 + exp x)) . VS.toList

    let showRes [dv,h1,h2] =
                 "D = "  ++ showFFloat (Just 3) dv ", " ++
                 "H1 = " ++ showFFloat (Just 3) h1 ", " ++
                 "H2 = " ++ showFFloat (Just 3) h2 ""
        showRes _ = error "Wtf? (1)"

    -- Confidence interval:  PCA on Hessian matrix, then for each
    -- eigenvalue λ add/subtract 1.96 / sqrt λ times the corresponding
    -- eigenvalue to the estimate.  Should describe a nice spheroid.
    let D2 _val grd hss = llk2 tab (paramVector2 $ VS.toList fit)
        d               = U.length grd
        (evals, evecs)  = eigSH' $ (d >< d) (U.toList hss)
        intervs         = [ (xform (fit + scale lam evec), xform (fit + scale (-lam) evec))
                          | (eval, evec) <- zip (VS.toList evals) (toRows evecs), let lam = 1.96 / sqrt eval ]

    when verbose $ putStrLn $ unlines $
            (:) (show res ++ ", " ++ show stats { finalValue = finalValue stats - llk_rr }) $
            (:) (showRes $ xform fit) $
            map (\(u,v) -> "[ " ++ showRes u ++ " .. " ++ showRes v ++ " ]") intervs

    return $! DivEst (xform fit) intervs

llk :: U.Vector Int -> [AD] -> AD
llk tab [delta,eta,eta2] = llk' tab 0 delta eta + llk' tab 6 delta eta2
llk _ _ = error "Wtf? (3)"

llk2 :: U.Vector Int -> [AD2] -> AD2
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

    !g_RR =        1 / Pr (log1p (exp delta))
    !g_AA = Pr delta / Pr (log1p (exp delta)) *      1 / Pr (log1p (exp eta))
    !g_RA = Pr delta / Pr (log1p (exp delta)) * Pr eta / Pr (log1p (exp eta))

    block ix g1 g2 g3 = U.ifoldl' step 0 $ U.slice (ix * maxD2) maxD2 tab
      where
        step !acc !i !num = acc - fromIntegral num * unPr p
          where
            (!d1,!d2) = i `quotRem` maxD
            p = g1 + Pr (- fromIntegral d1) * g2 + Pr (- fromIntegral (d1+d2)) * g3

