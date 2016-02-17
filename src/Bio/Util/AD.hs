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
    recip (C x)   = C (recip x)
    recip (D x u) = D (recip x) (U.map (y*) u) where y = negate $ recip $ x*x

    {-# INLINE fromRational #-}
    fromRational = C . fromRational


instance Floating AD where
    {-# INLINE pi #-}
    pi = C pi

    {-# INLINE exp #-}
    exp (C x)   = C (exp x)
    exp (D x u) = D (exp x) (U.map (* exp x) u)

    {-# INLINE sqrt #-}
    sqrt (C x)   = C (sqrt x)
    sqrt (D x u) = D (sqrt x) (U.map (* w) u) where w = recip $ 2 * sqrt x

    {-# INLINE log #-}
    log (C x)   = C (log x)
    log (D x u) = D (log x) (U.map (* recip x) u)

    {-# INLINE sin #-}
    sin (C x)   = C (sin x)
    sin (D x u) = D (sin x) (U.map (* cos x) u)

    {-# INLINE cos #-}
    cos (C x)   = C (cos x)
    cos (D x u) = D (cos x) (U.map (* negate (sin x)) u)

    {-
    tan = undefined -- :: a -> a
    asin = undefined -- :: a -> a
    atan = undefined -- :: a -> a
    acos = undefined -- :: a -> a
    sinh = undefined -- :: a -> a
    tanh = undefined -- :: a -> a
    cosh = undefined -- :: a -> a
    asinh = undefined -- :: a -> a
    atanh = undefined -- :: a -> a
    acosh = undefined -- :: a -> a -}


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
quietParameters = defaultParameters { printFinal = False, verbose = Quiet, maxItersFac = 20 }

debugParameters :: Parameters
debugParameters = defaultParameters { verbose = Verbose }

