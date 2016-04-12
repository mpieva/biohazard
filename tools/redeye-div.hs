-- Goal here is to read the tables from sample_div_tables, add them up, estimate divergence and
-- heterozygosity from them, store it back.

import Bio.Genocall.Metadata
import Bio.Util.AD
import Bio.Util.AD2
import Bio.Util.Numeric              ( log1p )
import Numeric                       ( showFFloat )
import Numeric.LinearAlgebra.HMatrix ( eigSH', (><), toRows, scale )

import qualified Data.Vector.Storable           as VS
import qualified Data.Vector.Unboxed            as U

-- XXX  If sample_divergences is set, we use it.  It it isn't, but we have sample_div_tables, we
-- estimate it, store it(!), then use it.  If we have neither, we bail.
{-
                conf_report $ "Estimating divergence parameters for " ++ sample ++ "..."
                est <- uncurry estimateSingle (foldr ... sample_div_tables)
                updateMetadata (H.adjust (\smp' -> smp' { sample_divergences = Just est })
                               (fromString sample)) (fromString conf_metadata)
-}


-- XXX need to think about what to return or store and why...
-- XXX we should estimate an indel rate, to be appended as the fourth result
estimateSingle :: Double -> U.Vector Int -> IO [Double]
estimateSingle llk_rr tab = do
    (fit, res, stats) <- minimize quietParameters 0.0001 (llk tab) (U.fromList [0,0,0])
    putStrLn $ show res ++ ", " ++ show stats { finalValue = finalValue stats - llk_rr }

    let showRes xx =
          case VS.toList xx of
            [delta, eta, eta2] ->
              let dv  = exp delta / (1 + exp delta)
                  ht  = exp eta  / (1 + exp eta) -- * dv
                  ht' = exp eta2 / (1 + exp eta2) -- * dv
              in "D = " ++ showFFloat (Just 3) dv ", " ++
                 "H1 = " ++ showFFloat (Just 3) ht ", " ++
                 "H2 = " ++ showFFloat (Just 3) ht' []
            _ -> error "Wtf? (1)"

    -- confidence interval:  PCA on Hessian matrix, then for each
    -- eigenvalue λ add/subtract 1.96 / sqrt λ times the corresponding
    -- eigenvalue to the estimate.  That should describe a nice
    -- spheroid.

    let D2 _val grd hss = llk2 tab (paramVector2 $ VS.toList fit)
        d               = U.length grd
        (evals, evecs)  = eigSH' $ (d >< d) (U.toList hss)

    putStrLn $ showRes fit
    sequence_ [ putStrLn $ "[ " ++ showRes (fit + scale lam evec)
                      ++ " .. " ++ showRes (fit + scale (-lam) evec) ++ " ]"
              | (eval, evec) <- zip (VS.toList evals) (toRows evecs)
              , let lam = 1.96 / sqrt eval ]

    case VS.toList fit of
        [delta, eta, eta2] ->
              let dv  = exp delta / (1 + exp delta)
                  ht  = exp eta  / (1 + exp eta) -- * dv
                  ht' = exp eta2 / (1 + exp eta2) -- * dv
              in return [ dv, ht * dv, ht' * dv ]
        _ -> error "Wtf? (2)"


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

