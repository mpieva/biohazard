{-# LANGUAGE RecordWildCards, NamedFieldPuns, ScopedTypeVariables, BangPatterns, TemplateHaskell, MultiParamTypeClasses, TypeFamilies  #-}
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
-- work.  On top of that, everything seems totally broken.

import Bio.Base
import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Bam.Rec
import Bio.Genocall.Adna
import Bio.Iteratee
import Control.Applicative
import Data.Bits
import Data.Foldable
import Data.Monoid
import Data.Traversable
import Data.Vec ( vec, Mat44, dot, getElem )
import Numeric.AD
import Numeric.AD.Internal.Forward.Double
import Numeric.Optimization.Algorithms.HagerZhang05
import Data.Time.Clock

import qualified Data.Vec as Vec
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Fusion.Stream as Stream
import qualified Numeric.AD.Mode.Forward.Double as Fwd

import Prelude hiding ( mapM_ )
import Data.Vector.Unboxed.Deriving

-- Log-likelihood of a bunch of data.  First argument is the bunch, second
-- is the set of parameters.
--
-- Let's see how this works out with the parameters.  Might want to
-- replace probabilities with log-likelihood-ratios or similar.

data Parms a = Parms !a !a !a !a !a !a

fromVec :: U.Vector Double -> Parms Double
fromVec vec = Parms (vec `U.unsafeIndex` 0) (vec `U.unsafeIndex` 1) (vec `U.unsafeIndex` 2)
                    (vec `U.unsafeIndex` 3) (vec `U.unsafeIndex` 4) (vec `U.unsafeIndex` 5)

toVec :: Parms Double -> U.Vector Double
toVec (Parms u v w x y z) = U.fromListN 6 [u,v,w,x,y,z]

{-# INLINE plus #-}
plus :: Num a => Parms a -> Parms a -> Parms a
Parms a b c d e f `plus` Parms u v w x y z = Parms (a+u) (b+v) (c+w) (d+x) (e+y) (f+z)

instance Functor Parms where
    fmap f (Parms u v w x y z) = Parms (f u) (f v) (f w) (f x) (f y) (f z)

instance Foldable Parms where
    foldr f n (Parms u v w x y z) = f u $! f v $! f w $! f x $! f y $! f z n
    foldl' f a (Parms u v w x y z) = (f $! (f $! (f $! (f $! (f $! f a u) v) w) x) y) z

instance Traversable Parms where
    traverse f (Parms u v w x y z) = Parms <$> f u <*> f v <*> f w <*> f x <*> f y <*> f z

sigmoid, isigmoid :: (Num a, Fractional a, Floating a) => a -> a
sigmoid x = 0.5 + 0.5 * x / sqrt (1+x*x)
isigmoid z = y / sqrt (1-y*y) where y = 2*z-1

{-# INLINE lk_fun1 #-}
lk_fun1 :: forall a . (Num a, Fractional a, Floating a) => U.Vector Word8 -> Parms a -> a
lk_fun1 bb parms = negate $ log $ lk bb
  where
    p_self :: a
    p_subst = (1 - p_self) * 0.333

    Parms l_endo l_subst l_sigma l_delta l_lam l_kap = parms

    f_endo = sigmoid l_endo
    p_self = sigmoid l_subst
    ssd_sigma = sigmoid l_sigma
    ssd_delta = sigmoid l_delta
    ssd_lambda = exp l_lam
    ssd_kappa = exp l_kap

    -- Technically, its transpose.  But it's symmetric anyway.
    subst_mat :: Mat44 a
    subst_mat = vec4 (vec4 p_self p_subst p_subst p_subst)
                     (vec4 p_subst p_self p_subst p_subst)
                     (vec4 p_subst p_subst p_self p_subst)
                     (vec4 p_subst p_subst p_subst p_self)

    lk :: U.Vector Word8 -> a
    lk br = {-f_endo *-} Stream.foldl' (*) 1 (Stream.zipWith lk1 (G.stream br) $
                                          G.stream (ssDamage SSD{..} False (U.length br)))
            -- + (1-f_endo) * G.foldl' (\acc pr -> acc * lk1 pr Vec.identity) 1 br

    lk1 :: Word8 -> Mat44 a -> a
    lk1 pr m | pr > 15   = 1
             | otherwise = getElem b m `dot` getElem a subst_mat
      where
        a = fromIntegral $ pr `shiftR` 2  -- from
        b = fromIntegral $ pr .&. 3       -- to


lkfun :: V.Vector (U.Vector Word8) -> U.Vector Double -> Double
lkfun brs parms = let !p' = fromVec parms
                  in V.foldl' (\a b -> a + lk_fun1 b p') 0 brs

gradfn :: V.Vector (U.Vector Word8) ->  U.Vector Double -> U.Vector Double
gradfn brs parms = let !p' = fromVec parms
                   in toVec $ V.foldl' (\a b -> a `plus` grad (lk_fun1 b) p')
                                       (Parms 0 0 0 0 0 0) brs

combofn :: V.Vector (U.Vector Word8) -> U.Vector Double -> (Double, U.Vector Double)
combofn brs parms = let !p' = fromVec parms
                        (!z,!zg) = V.foldl' (\(!a,!ag) bb -> let (!b,!bg) = grad' (lk_fun1 bb) p'
                                                             in (a + b, ag `plus` bg))
                                            (0, Parms 0 0 0 0 0 0) brs
                    in (z, toVec zg)

main = do
    brs <- concatDefaultInputs >=> run $ \_ ->
           joinI $ filterStream (not . br_isUnmapped) $
           joinI $ mapStream pack_record $
           joinI $ filterStream (U.all (<16)) $
           stream2vectorN 1000


    let v0 = U.fromList [3,(-5),2,(-4),0,0]

    print $ V.length brs
    print' v0
    print =<< getCurrentTime
    print $ lkfun brs v0
    print =<< getCurrentTime
    print $ gradfn brs v0
    print =<< getCurrentTime
    print $ combofn brs v0
    print =<< getCurrentTime

    let params = defaultParameters { verbose = VeryVerbose }

    (xs, r, st) <- optimize params 0.0001 v0
                            (VFunction $ lkfun brs)
                            (VGradient $ gradfn brs)
                            (Just . VCombined $ combofn brs)

    print' xs
    print r
    print st


print' vec =
    print $ map sigmoid (G.toList $ G.take 4 vec)
         ++ map exp (G.toList $ G.drop 4 vec)


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
                               2 -> mk_pair' 4
                               4 -> mk_pair' 8
                               8 -> mk_pair' 12
                               _ -> const esc

    mk_pair' :: Word8 -> Nucleotides -> Word8
    mk_pair' a (Ns b) = case b of 1 -> a .|. 0
                                  2 -> a .|. 1
                                  4 -> a .|. 2
                                  8 -> a .|. 3
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




