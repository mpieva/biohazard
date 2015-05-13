{-# LANGUAGE OverloadedStrings, TemplateHaskell, BangPatterns, RecordWildCards, FlexibleContexts, TypeFamilies #-}
-- Scan file with GT likelihoods, fit something...
--
-- First iteration:  Mitochondrion only.   We don't need to fit
-- anything.  So far, the likelihoods behave strangely in that smaller
-- \theta is always better, as long as it doesn't become zero.
--
-- Second iteration:  Mitochondrion only, but with a divergence
-- parameter.  Needs to be scanned in parallel with a TwoBit file.

import Bio.Bam.Pileup
import Bio.Genocall.AvroFile
import Bio.Iteratee
import Bio.TwoBit
import Bio.Util
import Control.Monad ( zipWithM_, forM_, foldM )
import Control.Monad.Primitive
import Data.Avro
import Data.ByteString ( ByteString )
import Data.HashMap.Strict ( toList )
import Data.List ( foldl' )
import Data.MiniFloat ( mini2float )
import Data.Text.Encoding ( encodeUtf8, decodeUtf8 )
import Data.Text ( Text, unpack )
import Data.Strict.Tuple
import Numeric.Optimization.Algorithms.HagerZhang05

import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM

import AD

{-
main = do -- let ctr :: ByteString
          ctr <- enumPure1Chunk test_data >=> run $ joinI $ writeAvroContainer ctr_opts $ stream2stream
          print ctr
          rd <- enumPure1Chunk (ctr :: ByteString) >=> run $ joinI $
                (readAvroContainer {- :: Enumeratee ByteString [Int] IO [Int] -}) $ stream2list

          print rd >> print (rd == test_data)

  where
    ctr_opts = ContainerOpts 1 "test"

    test_data :: [GenoCallBlock]
    test_data = [ GenoCallBlock "chr1" 0
                    [ GenoCallSite (CallStats 99 50 3500 99000)
                                   [0,1,2]
                                   (CallStats 80 40 3000 47000)
                                   [ IndelVariant 1 (V_Nuc $ fromList $ read "ACGT") ]
                                   [1001,1002,1003] ]
                , GenoCallBlock "MT" 0 [] ]
-}

main :: IO ()
main = do hg19 <- openTwoBit "/mnt/datengrab/hg19.2bit"
          mtbl <- UM.replicate (max_lk-min_lk+1) 0
          p0 :!: pe <- enumDefaultInputs >=> run $
                    joinI $ readAvroContainer $ \meta -> do
                        liftIO . forM_ (toList meta) $ \(k,v) ->
                            putStrLn $ unpack k ++ ": " ++ unpack (decodeUtf8 v)
                        foldStreamM (lk_block (all_lk mtbl) hg19) (1 :!: 1)

          tbl  <- U.unsafeFreeze mtbl
          print (unPr p0, p0)
          print (unPr pe, pe)
          mapM_ print [ U.slice i 64 tbl | i <- [ 0, 64 .. U.length tbl-1 ] ]
          print . unPr $ llk1 tbl pe (D 0.001 (U.fromList [1]))
  where
    all_lk tbl (p1 :!: p2) ref site = (lk0 p1 site :!:) `fmap` lk1 tbl p2 ref site

-- | Scans block together with reference sequence.  Folds a monadic
-- action over the called sites.
lk_block :: Monad m => (b -> Nucleotide -> GenoCallSite -> m b) -> TwoBitFile -> b -> GenoCallBlock -> m b
lk_block f tbf b GenoCallBlock{..} = foldM2 b refseq called_sites
  where
    refseq = getLazySubseq tbf (Pos (encodeUtf8 reference_name) start_position)

    foldM2 acc (x:xs) (y:ys) = do !acc' <- f acc x y ; foldM2 acc' xs ys
    foldM2 acc [    ]      _ = return acc
    foldM2 acc      _ [    ] = return acc


-- | Likelihood with flat prior (no parameters).
lk0 :: Prob Double -> GenoCallSite -> Prob Double
lk0 !pp GenoCallSite{..} | U.length snp_likelihoods == 4 =
    pp * 0.25 * U.sum (U.map (Pr . negate . mini2float) snp_likelihoods)

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
type LkTable  = U.Vector  Int

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

-- | Actual log-likelihood.  Gets a table and a divergence value.
-- Returns likelihoods and first two derivatives with respect to the
-- divergence value.
-- llk :: (Ord a, Floating a) => LkTable -> a -> Prob a
llk1 :: LkTable -> Prob Double -> AD -> Prob AD
llk1 tbl (Pr p) d = Pr (C p) * U.ifoldl' step 1 tbl
  where
    !d1 = toProb $ 1-d
    !d3 = toProb $ d/3

    step acc qq num = acc * ( d1 + d3 * Pr (fromIntegral $ min_lk + qq) ) `pow` num

