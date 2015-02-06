{-# LANGUAGE RecordWildCards, NamedFieldPuns, BangPatterns, MultiParamTypeClasses, TypeFamilies  #-}
-- Estimates aDNA damage.  Crude first version.
--
-- - Read a BAM file, make compact representation of the reads.
-- - Compute likelihood of read under simple model of
--   damage, error/divergence, contamination.
-- - Iterate using Data.AD.Newton or similar.
--
-- For the less crude version:  follow the iobio subsampling strategy
-- using an index.

-- Trying to compute symbolically is too much, the high power terms get
-- out of hand quickly, and we get mixed powers of \lambda and \kappa.

-- For the fitting, we simplify radically: ignore sequencing error,
-- assume damage and simple, symmetric substitutions.

-- So far, the fastest version uses a custom data type (Parms), which we
-- hope gets unboxed.  Unboxed Vectors don't work, because they aren't
-- polymorphic enough.  The optimizer is something external.
-- Implementing AD ourselves (only the gradient, forward mode) could
-- work.
--
-- If I include parameters, whose true value is zero, the transformation
-- to the log-odds-ratio doesn't work, because then the maximum doesn't
-- exist anymore.  A different transformation ('sigmoid2'/'isigmoid2'
-- below) allows for an actual zero (but not one), while avoiding ugly
-- boundary conditions.  That works.
--
-- The current hack assumes all molecules have an overhang at both ends,
-- then each base gets deaminated with a position dependent probability.
-- If we try to model a fraction of undeaminated molecules in addition,
-- this fails.  To rescue the idea, I guess we must really decide if the
-- molecule has an overhang at all (probability 1/2) at each end, then
-- deaminate it.

-- TODO:
-- Before this can be packaged and used, the following needs to be done:
-- - Subsample large BAM files (don't always take the beginning),
-- - Start with a crude estimate of parameters,
-- - Fix the model(s), so a contaminant fraction can be estimated.
-- - Implement both SSD and DSD.


import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Bam.Rec
import Bio.Base
import Bio.Genocall.Adna
import Bio.Iteratee
import Control.Applicative
import Data.Bits
import Data.Foldable
import Data.Monoid
import Data.Traversable
import Data.Vec ( vec, Mat44, dot, getElem )
import Numeric.Optimization.Algorithms.HagerZhang05

import qualified Data.Vector                as V
import qualified Data.Vector.Fusion.Stream  as S
import qualified Data.Vector.Generic        as G
import qualified Data.Vector.Unboxed        as U

import AD
import Prelude hiding ( mapM_ )

sigmoid, sigmoid2, isigmoid, isigmoid2 :: (Num a, Fractional a, Floating a) => a -> a
sigmoid l = 1 / ( 1 + exp l )
isigmoid p = log $ ( 1 - p ) / p

sigmoid2 x = y*y where y = (exp x - 1) / (exp x + 1)
isigmoid2 y = log $ (1 + sqrt y) / (1 - sqrt y)

{-# INLINE lk_fun1 #-}
lk_fun1 :: (Num a, Fractional a, Floating a) => U.Vector Word8 -> [a] -> a
lk_fun1 bb parms = negate $ log $ lk bb
  where
    l_subst:l_sigma:l_delta:l_lam:l_kap:_ = parms

    -- Good initial guesses may be necessary, too.
    -- f_exo = sigmoid l_endo
    p_self = 1 - sigmoid2 l_subst
    p_subst = (1 - p_self) * 0.333

    ssd_sigma  = sigmoid l_sigma
    ssd_delta  = sigmoid l_delta
    ssd_lambda = sigmoid l_lam
    ssd_kappa  = sigmoid l_kap

    -- Technically, its transpose.  But it's symmetric anyway.
    subst_mat = vec4 (vec4 p_self p_subst p_subst p_subst)
                     (vec4 p_subst p_self p_subst p_subst)
                     (vec4 p_subst p_subst p_self p_subst)
                     (vec4 p_subst p_subst p_subst p_self)

    lk br = {-(1-f_exo) *-} S.foldl' (*) 1 (S.zipWith lk1 (G.stream br) $
                                          G.stream (ssDamage SSD{..} False (U.length br)))
            -- + f_exo * G.foldl' (\acc pr -> acc * lk1 pr W.identity) 1 br

    lk1 pr m | pr > 15   = 1
             | otherwise = getElem b m `dot` getElem a subst_mat
      where
        a = fromIntegral $ pr `shiftR` 2  -- from
        b = fromIntegral $ pr .&. 3       -- to


lkfun :: V.Vector (U.Vector Word8) -> U.Vector Double -> Double
lkfun brs parms = V.foldl' (\a b -> a + lk_fun1 b ps) 0 brs
  where
    !ps = U.toList parms

combofn :: V.Vector (U.Vector Word8) -> U.Vector Double -> (Double, U.Vector Double)
combofn brs parms = (x,g)
  where
    !ps     = paramVector $ U.toList parms
    (D x g) = V.foldl' (\a b -> a + lk_fun1 b ps) 0 brs


main :: IO ()
main = do
    brs <- concatDefaultInputs >=> run $ \_ ->
           joinI $ filterStream (not . br_isUnmapped) $
           joinI $ mapStream pack_record $
           joinI $ filterStream (U.all (<16)) $
           stream2vectorN 1000

    mapM_ print brs

    let v0 = U.fromList $ {-isigmoid2 0.05 :-} isigmoid2 0.001 : map isigmoid [0.02, 0.5, 0.3, 0.3]

    print $ V.length brs
    print' v0
    print $ lkfun brs v0
    print $ snd $ combofn brs v0
    print $ combofn brs v0

    let params = defaultParameters { verbose = VeryVerbose }

    (xs, r, st) <- optimize params 0.0000000001 v0
                            (VFunction $ lkfun brs)
                            (VGradient $ snd . combofn brs)
                            (Just . VCombined $ combofn brs)

    print' xs
    print r
    print st


print' vec = print $
    map sigmoid2 (G.toList $ G.take 1 vec) ++
    map sigmoid (G.toList $ G.drop 1 vec)


-- We'll require the MD field to be present.  Then we cook each read
-- into a list of paired bases.  Deleted bases are dropped, inserted
-- bases replaced with an escape code.
--
-- XXX  This is annoying... almost, but not quite the same as the code
-- in the "Pileup" module.  This also relies on MD and doesn't offer the
-- alternative of accessing a reference genome.  The latter may not be
-- worth the trouble.

pack_record :: BamRaw -> U.Vector Word8
pack_record br = if br_isReversed br then revcom u1 else u1
  where
    BamRec{..} = decodeBamEntry br

    revcom = U.reverse . U.map (xor 15)

    u1 = U.fromList $ go (unCigar b_cigar) (U.toList b_seq) (br_get_md br)

    esc = 16

    mk_pair :: Nucleotides -> Nucleotides -> Word8
    mk_pair (Ns a) = case a of 1 -> mk_pair' 0
                               2 -> mk_pair' 1
                               4 -> mk_pair' 2
                               8 -> mk_pair' 3
                               _ -> const esc

    mk_pair' :: Word8 -> Nucleotides -> Word8
    mk_pair' a (Ns b) = case b of 1 -> a .|. 0
                                  2 -> a .|. 4
                                  4 -> a .|. 8
                                  8 -> a .|. 12
                                  _ -> esc

    go :: [(CigOp,Int)] -> [Nucleotides] -> [MdOp] -> [Word8]

    go ((_,0):cs)   ns mds  = go cs ns mds
    go cs ns (MdNum  0:mds) = go cs ns mds
    go cs ns (MdDel []:mds) = go cs ns mds
    go  _ []              _ = []

    go ((Mat,nm):cs) (n:ns) (MdNum mm:mds) = mk_pair n n  : go ((Mat,nm-1):cs) ns (MdNum (mm-1):mds)
    go ((Mat,nm):cs) (n:ns) (MdRep n':mds) = mk_pair n n' : go ((Mat,nm-1):cs) ns               mds
    go ((Mat,nm):cs)    ns  (MdDel ds:mds) =                go ((Mat, nm ):cs) ns               mds

    go ((Ins,nm):cs) ns mds = replicate nm esc ++ go cs (drop nm ns) mds
    go ((SMa,nm):cs) ns mds = replicate nm esc ++ go cs (drop nm ns) mds
    go ((Del,nm):cs) ns (MdDel (_:ds):mds) = go ((Del,nm-1):cs) ns (MdDel ds:mds)
    go ((Del,nm):cs) ns (           _:mds) = go ((Del, nm ):cs) ns           mds

    go (_:cs) nd mds = go cs nd mds




