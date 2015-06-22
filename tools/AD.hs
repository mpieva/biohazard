{-# LANGUAGE BangPatterns #-}
module AD ( AD(..), liftD, sqr, lsum, llerp
          , paramVector, minimize
          , module Numeric.Optimization.Algorithms.HagerZhang05
          ) where

import Bio.Util ( log1p )
import Data.List ( foldl1' )
import Numeric.Optimization.Algorithms.HagerZhang05
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Storable as V

-- Simple forward-mode AD to get a scalar valued function and a
-- gradient.

data AD = C !Double | D !Double !(U.Vector Double)
  deriving Show

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
    sqrt (D x u) = D (sqrt x) (U.map (*w) u) where w = recip $ 2 * sqrt x

    {-# INLINE log #-}
    log (C x)   = C (log x)
    log (D x u) = D (log x) (U.map (*w) u) where w = recip x

    {- (**) = undefined -- :: a -> a -> a
    logBase = undefined -- :: a -> a -> a
    sin = undefined -- :: a -> a
    tan = undefined -- :: a -> a
    cos = undefined -- :: a -> a
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

-- | Computes \( \log ( \sum_i e^{x_i} ) \) sensibly.  The list must be
-- sorted in descending(!) order.
{-# INLINE lsum #-}
lsum :: [AD] -> AD
lsum xs = foldl1' (\x y -> if x >= y then x + log1p (exp (y-x)) else err) xs
    where err = error $ "lsum: argument list must be in descending order: " ++ show xs

-- | Computes \( \log \left( c e^x + (1-c) e^y \right) \).
{-# INLINE llerp #-}
llerp :: AD -> AD -> AD -> AD
llerp c x y | c == 0.0  = y
            | c == 1.0  = x
            | x >= y    = log     c  + x + log1p ( (1-c)/c * exp (y-x) )
            | otherwise = log1p (-c) + y + log1p ( c/(1-c) * exp (x-y) )

