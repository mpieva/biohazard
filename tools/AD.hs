{-# LANGUAGE BangPatterns #-}
module AD where

import qualified Data.Vector.Unboxed as U

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


{-
-- | Number and first two derivatives.  Only one argument for the time
-- being.
data AD3 = AD3 !Double !Double !Double

instance Num AD3 where
    {-# INLINE (+) #-}
    AD3 x u a + AD3 y v b = AD3 (x+y) (u+v) (a+b)

    {-# INLINE (-) #-}
    AD3 x u a - AD3 y v b = AD3 (x-y) (u-v) (a-b)

    {-# INLINE (*) #-}
    AD3 x u a * AD3 y v b = AD3 (x*y) (u*y + x*v) (a*v + u*b)

    {-# INLINE negate #-}
    negate (AD3 y v b) = AD3 (negate y) (negate v) (negate b)

    {-# INLINE fromInteger #-}
    fromInteger x = AD3 (fromInteger x) 0 0

    {-# INLINE abs #-}
    abs (AD3 x u a) | x < 0     = AD3 (negate x) (negate u) (negate a)
                    | otherwise = AD3 x u a

    {-# INLINE signum #-}
    signum (AD3 x _ _)   = AD3 (signum x) 0 0

instance Fractional AD3 where
    {-# INLINE (/) #-}
    AD3 x u a / AD3 y v b = AD3 (x/y) ((u*y-x*v)/(y*y)) ((a*v-u*b)/(v*v))

    {-# INLINE recip #-}
    recip (AD3 x u a) = AD3 (recip x) (-u/(x*x)) (-a/(u*u))

    {-# INLINE fromRational #-}
    fromRational x = AD3 (fromRational x) 0 0


instance Floating AD where
    {-# INLINE pi #-}
    pi = AD3 pi 0 0

    {-# INLINE exp #-}
    exp (AD3 x u a) = AD3 (exp x) (exp x * u) (
    exp (C x)   = C (exp x)
    exp (D x u) = D (exp x) (U.map (* exp x) u)

    {-# INLINE sqrt #-}
    sqrt (C x)   = C (sqrt x)
    sqrt (D x u) = D (sqrt x) (U.map (*w) u) where w = recip $ 2 * sqrt x

    {-# INLINE log #-}
    log (C x)   = C (log x)
    log (D x u) = D (log x) (U.map (*w) u) where w = recip x
-}

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

