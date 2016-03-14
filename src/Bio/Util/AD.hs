{-# LANGUAGE BangPatterns #-}
module Bio.Util.AD
          ( AD(..), paramVector, minimize
          , module Numeric.Optimization.Algorithms.HagerZhang05
          , debugParameters, quietParameters
          ) where

import Numeric.Optimization.Algorithms.HagerZhang05
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Storable as V

-- | Simple forward-mode AD to get a scalar valued function with gradient.
data AD = C !Double | D !Double !(U.Vector Double) deriving Show

instance Eq AD where
    C x   == C y   = x == y
    C x   == D y _ = x == y
    D x _ == C y   = x == y
    D x _ == D y _ = x == y

instance Ord AD where
    C x   `compare` C y   = x `compare` y
    C x   `compare` D y _ = x `compare` y
    D x _ `compare` C y   = x `compare` y
    D x _ `compare` D y _ = x `compare` y

instance Num AD where
    {-# INLINE (+) #-}
    C x   + C y   = C (x+y)
    C x   + D y v = D (x+y) v
    D x u + C y   = D (x+y) u
    D x u + D y v = D (x+y) (U.zipWith (+) u v)

    {-# INLINE (-) #-}
    C x   - C y   = C (x-y)
    C x   - D y v = D (x-y) (U.map negate v)
    D x u - C y   = D (x-y) u
    D x u - D y v = D (x-y) (U.zipWith (-) u v)

    {-# INLINE (*) #-}
    C x   * C y   = C (x*y)
    C x   * D y v = D (x*y) (U.map (x*) v)
    D x u * C y   = D (x*y) (U.map (y*) u)
    D x u * D y v = D (x*y) (U.zipWith (+) (U.map (x*) v) (U.map (y*) u))

    {-# INLINE negate #-}
    negate (C x)   = C (negate x)
    negate (D x u) = D (negate x) (U.map negate u)

    {-# INLINE fromInteger #-}
    fromInteger = C . fromInteger

    {-# INLINE abs #-}
    abs (C x) = C (abs x)
    abs (D x u) | x < 0     = D (negate x) (U.map negate u)
                | otherwise = D x u

    {-# INLINE signum #-}
    signum (C x)   = C (signum x)
    signum (D x _) = C (signum x)


instance Fractional AD where
    {-# INLINE (/) #-}
    C x   / C y   = C (x/y)
    D x u / C y   = D (x*z) (U.map (z*) u) where z = recip y
    C x   / D y v = D (x/y) (U.map (w*) v) where w = negate $ x * z * z ; z = recip y
    D x u / D y v = D (x/y) (U.zipWith (-) (U.map (z*) u) (U.map (w*) v))
        where z = recip y ; w = x * z * z

    {-# INLINE recip #-}
    recip = liftF recip (\x -> - recip (x*x))

    {-# INLINE fromRational #-}
    fromRational = C . fromRational


instance Floating AD where
    {-# INLINE pi #-}
    pi = C pi

    {-# INLINE exp #-}
    exp   = liftF exp exp

    {-# INLINE sqrt #-}
    sqrt  = liftF sqrt $ \x -> recip (2 * sqrt x)

    {-# INLINE log #-}
    log   = liftF log recip

    sin   = liftF sin cos
    cos   = liftF cos (negate . sin)
    sinh  = liftF sinh cosh
    cosh  = liftF cosh sinh

    tan   = liftF tan   $ \x ->   recip (cos x * cos x)
    tanh  = liftF tanh  $ \x ->   recip (cosh x * cosh x)
    asin  = liftF asin  $ \x ->   recip (sqrt (1 - x * x))
    acos  = liftF acos  $ \x -> - recip (sqrt (1 - x * x))
    atan  = liftF atan  $ \x ->   recip (1 + x * x)
    asinh = liftF asinh $ \x ->   recip (sqrt (x * x + 1))
    acosh = liftF acosh $ \x -> - recip (sqrt (x * x - 1))
    atanh = liftF atanh $ \x ->   recip (1 - x * x)


{-# INLINE liftF #-}
liftF :: (Double -> Double) -> (Double -> Double) -> AD -> AD
liftF f _ (C x) = C (f x)
liftF f g (D x u) = D (f x) (U.map (* g x) u)

{-# INLINE paramVector #-}
paramVector :: [Double] -> [AD]
paramVector xs = [ D x (U.generate l (\j -> if i == j then 1 else 0)) | (i,x) <- zip [0..] xs ]
  where l = length xs

{-# INLINE minimize #-}
minimize :: Parameters -> Double -> ([AD] -> AD) -> U.Vector Double -> IO (V.Vector Double, Result, Statistics)
minimize params eps func v0 =
    optimize params eps v0 (VFunction  $ fst . combofn)
                           (VGradient  $ snd . combofn)
                           (Just . VCombined $ combofn)
  where
    combofn parms = case func $ paramVector $ U.toList parms of
                D x g -> ( x, g )
                C x   -> ( x, U.replicate (U.length parms) 0 )


quietParameters :: Parameters
quietParameters = defaultParameters { printFinal = False, verbose = Quiet, maxItersFac = 123 }

debugParameters :: Parameters
debugParameters = defaultParameters { verbose = Verbose }

