{-# LANGUAGE OverloadedStrings, TemplateHaskell, BangPatterns, RecordWildCards, FlexibleContexts, TypeFamilies #-}
-- Scan file with GT likelihoods, fit something...
--
-- First iteration:  Mitochondrion only.   We don't need to fit
-- anything.  So far, the likelihoods behave strangely in that smaller
-- \theta is always better, as long as it doesn't become zero.
--
-- Second iteration:  Mitochondrion only, but with a divergence
-- parameter.  Needs to be scanned in parallel with a TwoBit file.

import Bio.Genocall.AvroFile
import Bio.Iteratee
import Bio.TwoBit
import Bio.Util
import Control.Monad.Primitive
import Data.Avro
import Data.List ( intercalate )
import Data.MiniFloat ( mini2float )
import Data.Text.Encoding ( encodeUtf8 )
import Data.Strict.Tuple ( Pair((:!:)), (:!:) )
import Numeric ( showFFloat )
import Numeric.Optimization.Algorithms.HagerZhang05

import qualified Data.Vector.Storable as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM

import AD

main :: IO ()
main = do hg19 <- openTwoBit "/mnt/datengrab/hg19.2bit"
          mtbl <- UM.replicate (max_lk-min_lk+1) 0

          let all_lk tbl (p1 :!: p2) ref site = (lk0 p1 site :!:) `fmap` lk1 tbl p2 ref site

          p0 :!: pe <- enumDefaultInputs >=> run $
                    joinI $ readAvroContainer $ \_meta -> do
                        -- liftIO . forM_ (toList meta) $ \(k,v) ->
                            -- putStrLn $ unpack k ++ ": " ++ unpack (decodeUtf8 v)
                        foldStreamM (lk_block (all_lk mtbl) hg19) (1 :!: 1)

          tbl  <- U.unsafeFreeze mtbl

          -- optimize llk1 vs. d argument.
          let plainfn :: U.Vector Double -> Double
              plainfn args = llk1 tbl (unPr pe) $ args U.! 0

              combofn :: U.Vector Double -> (Double, U.Vector Double)
              combofn args = case llk1 tbl (C (unPr pe)) $ D (args U.! 0) (U.singleton 1) of
                                (D x dx) -> ( x, dx )

              params = defaultParameters { printFinal = False, verbose = {- Verbose -} Quiet, maxItersFac = 20 }

          (x,q,s) <- optimize params 0.0001 (U.singleton 0.01)
                            (VFunction plainfn)
                            (VGradient $ snd . combofn)
                            (Just $ VCombined combofn)

          -- print $ llk1 tbl (C (unPr pe)) (D 0.001 (U.singleton 1))
          putStrLn $ intercalate "\t"
            [ showNum . round $ unPr p0, showFFloat (Just 5) (sigmoid2 $ x V.! 0) []
            , show q, showNum . round $ finalValue s, show s ]
          -- print (map sigmoid2 $ V.toList x, q, s)


-- | Scans block together with reference sequence.  Folds a monadic
-- action over the called sites.
lk_block :: Monad m => (b -> Nucleotide -> GenoCallSite -> m b) -> TwoBitFile -> b -> GenoCallBlock -> m b
lk_block f tbf b GenoCallBlock{..} = foldM3f b start_position refseq called_sites
  where
    refseq = getLazySubseq tbf (Pos (encodeUtf8 reference_name) start_position)

    foldM2 acc (x:xs) (y:ys) = do !acc' <- f acc x y ; foldM2 acc' xs ys
    foldM2 acc [    ]      _ = return acc
    foldM2 acc      _ [    ] = return acc


    -- XXX terrible hack to deal with PhiX!  Remove this as soon as
    -- sensible!
    foldM3f acc n (x:xs) (y:ys)
        | n `elem` bad          = foldM3f acc (succ n) xs ys
        | otherwise             = do !acc' <- f acc x y ; foldM3f acc' (succ n) xs ys
    foldM3f acc _ [    ]      _ = return acc
    foldM3f acc _      _ [    ] = return acc

    bad = [1400,1643]
    -- bad = [586,832,1649,2810,4517]

{- p_block tbf GenoCallBlock{..} = do
    printf "Block %s:%d-%d\n" (show reference_name) start_position
                              (start_position + length called_sites)
    zipWithM_ (curry print) refseq (map snp_likelihoods called_sites)
  where
    refseq = getLazySubseq tbf (Pos (encodeUtf8 reference_name) start_position) -}



-- | Likelihood with flat prior (no parameters).
lk0 :: Prob Double -> GenoCallSite -> Prob Double
lk0 !pp GenoCallSite{..} | U.length snp_likelihoods == 4 =
    pp * 0.25 * U.sum (U.map (Pr . negate . mini2float) snp_likelihoods)
                         | otherwise = pp

-- | Likelihood precomputation.  Total likelihood computes as product
-- over sites @i@ with reference alleles @X_i@:
-- @
--   L(d) = \prod_i ( (1-d) * GL(X_i) + 1/3 * d * \sum_{Y/=X_i} GL(Y) )
--        = \prod_i GL(X_i) * \prod_i ( 1 - d + 1/3 * d * \sum_{Y/=X_i} GL(Y)/GL(X) )
-- @
--
-- We compute the first term on the first pass and tabulate a quantized
-- form of the second term: @round (log \sum_{Y/=X_i} GL(Y)/GL(X))@.
-- (Maybe add a scaling factor, though the plain natural log seems
-- pretty good.)

type LkTableM = UM.MVector (PrimState IO) Int
type LkTable  = U.Vector                  Int

min_lk, max_lk :: Int
min_lk = -256
max_lk =  255

-- | Likelihood with one parameter, the divergence.  Computes one
-- part directly, bins the variable part into a mutable table.
lk1 :: LkTableM -> Prob Double -> Nucleotide -> GenoCallSite -> IO (Prob Double)
lk1 tbl !pp ref GenoCallSite{..} | U.length snp_likelihoods == 4 = do
    let lx   = Pr . negate . mini2float $ snp_likelihoods U.! fromEnum ref
        odds = U.ifoldl' (\a i v -> if i == fromEnum ref then a else a + Pr (- mini2float v)) 0 snp_likelihoods / lx
        qq   = round (unPr odds) `min` max_lk `max` min_lk   - min_lk
    UM.write tbl qq . succ =<< UM.read tbl qq
    return $! pp * lx
                                 | otherwise = return pp

-- | Actual negative log-likelihood.  Gets a table and a divergence
-- value.  Returns likelihoods and first two derivatives with respect to
-- the divergence value.
llk1 :: (Ord a, Floating a) => LkTable -> a -> a -> a
llk1 tbl p d = U.ifoldl' step (-p) tbl
  where
    !d1 = log1p (- sigmoid2 d)
    !d3 = log (sigmoid2 d) - log 3

    step acc qq num = acc - fromIntegral num * (d1 <#> d3 + fromIntegral (min_lk + qq))

