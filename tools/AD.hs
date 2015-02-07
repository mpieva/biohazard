{-# LANGUAGE BangPatterns #-}
module AD where

import Control.Monad.ST

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as M

-- Simple forward-mode AD to get a scalar valued function and a
-- gradient.

data AD = C !Double | D !Double !(U.Vector Double)
  deriving Show

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


paramVector :: [Double] -> [AD]
paramVector xs = [ D x (U.generate l (\j -> if i == j then 1 else 0)) | (i,x) <- zip [0..] xs ]
  where l = length xs

