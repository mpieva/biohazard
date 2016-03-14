{-# LANGUAGE BangPatterns #-}
module Bio.Util.AD2 ( AD2(..), paramVector2 ) where

import qualified Data.Vector.Unboxed as U

-- | Simple forward-mode AD to get a scalar valued function
-- with gradient and Hessian.
data AD2 = C2 !Double | D2 !Double !(U.Vector Double) !(U.Vector Double)

instance Show AD2 where
    show (C2 x) = show x
    show (D2 x y z) = show x ++ " " ++ show (U.toList y) ++ " "
                    ++ show [ U.toList (U.slice i d z) | i <- [0, d .. d*d-1] ]
        where d = U.length y

instance Eq AD2 where
    C2 x     == C2 y     = x == y
    C2 x     == D2 y _ _ = x == y
    D2 x _ _ == C2 y     = x == y
    D2 x _ _ == D2 y _ _ = x == y

instance Ord AD2 where
    C2 x     `compare` C2 y     = x `compare` y
    C2 x     `compare` D2 y _ _ = x `compare` y
    D2 x _ _ `compare` C2 y     = x `compare` y
    D2 x _ _ `compare` D2 y _ _ = x `compare` y

instance Num AD2 where
    {-# INLINE (+) #-}
    C2 x     + C2 y     = C2 (x+y)
    C2 x     + D2 y v h = D2 (x+y) v h
    D2 x u g + C2 y     = D2 (x+y) u g
    D2 x u g + D2 y v h = D2 (x+y) (U.zipWith (+) u v) (U.zipWith (+) g h)

    {-# INLINE (-) #-}
    C2 x     - C2 y     = C2 (x-y)
    C2 x     - D2 y v h = D2 (x-y) (U.map negate v) (U.map negate h)
    D2 x u g - C2 y     = D2 (x-y) u g
    D2 x u g - D2 y v h = D2 (x-y) (U.zipWith (-) u v) (U.zipWith (-) g h)

    {-# INLINE (*) #-}
    C2 x     * C2 y     = C2 (x*y)
    C2 x     * D2 y v h = D2 (x*y) (U.map (x*) v) (U.map (x*) h)
    D2 x u g * C2 y     = D2 (x*y) (U.map (y*) u) (U.map (y*) g)
    D2 x u g * D2 y v h = D2 (x*y) grad hess
      where grad = U.zipWith (+) (U.map (x*) v) (U.map (y*) u)
            hess = U.zipWith (+)
                        (U.zipWith (+) (U.map (x*) h) (U.map (y*) g))
                        (U.zipWith (+) (cross u v) (cross v u))

    {-# INLINE negate #-}
    negate (C2 x)     = C2 (negate x)
    negate (D2 x u g) = D2 (negate x) (U.map negate u) (U.map negate g)

    {-# INLINE fromInteger #-}
    fromInteger = C2 . fromInteger

    {-# INLINE abs #-}
    abs (C2 x) = C2 (abs x)
    abs (D2 x u g) | x < 0     = D2 (negate x) (U.map negate u) (U.map negate g)
                   | otherwise = D2 x u g

    {-# INLINE signum #-}
    signum (C2 x)     = C2 (signum x)
    signum (D2 x _ _) = C2 (signum x)


instance Fractional AD2 where
    {-# INLINE (/) #-}
    C2 x     / C2 y     = C2 (x/y)
    D2 x u g / C2 y     = D2 (x*z) (U.map (z*) u) (U.map (z*) g) where z = recip y
    x / y = x * recip y

    {-# INLINE recip #-}
    recip = liftF recip (\x -> - recip (sqr x)) (\x -> 2 * recip (cube x))

    {-# INLINE fromRational #-}
    fromRational = C2 . fromRational

instance Floating AD2 where
    {-# INLINE pi #-}
    pi = C2 pi

    {-# INLINE exp #-}
    exp = liftF exp exp exp

    {-# INLINE sqrt #-}
    sqrt = liftF sqrt (\x -> recip $ 2 * sqrt x) (\x -> - recip (sqrt (cube x)))

    {-# INLINE log #-}
    log = liftF log recip (\x -> - recip (sqr x))

    sin   = liftF sin cos (negate . sin)
    cos   = liftF cos (negate . sin) (negate . cos)
    sinh  = liftF sinh cosh sinh
    cosh  = liftF cosh sinh cosh

    tan   = liftF tan   (\x ->   recip (sqr (cos  x))) (\x ->  2 * tan  x / sqr (cos  x))
    tanh  = liftF tanh  (\x ->   recip (sqr (cosh x))) (\x -> -2 * tanh x / sqr (cosh x))
    
    asin  = liftF asin  (\x ->   recip (sqrt (1 - sqr x))) (\x ->      x / sqrt (cube (1 - sqr x)))
    acos  = liftF acos  (\x -> - recip (sqrt (1 - sqr x))) (\x ->     -x / sqrt (cube (1 - sqr x)))
    asinh = liftF asinh (\x ->   recip (sqrt (sqr x + 1))) (\x ->     -x / sqrt (cube (sqr x + 1)))
    acosh = liftF acosh (\x -> - recip (sqrt (sqr x - 1))) (\x ->      x / sqrt (cube (sqr x - 1)))
    atan  = liftF atan  (\x ->   recip       (1 + sqr x))  (\x -> -2 * x / sqr (1 + sqr x))
    atanh = liftF atanh (\x ->   recip       (1 - sqr x))  (\x ->  2 * x / sqr (1 - sqr x))

{-# INLINE sqr #-}
sqr :: Double -> Double
sqr x = x * x

{-# INLINE cube #-}
cube :: Double -> Double
cube x = x * x * x

{-# INLINE liftF #-}
liftF :: (Double -> Double) -> (Double -> Double) -> (Double -> Double) -> AD2 -> AD2
liftF f  _  _  (C2 x)     = C2 (f x)
liftF f f' f'' (D2 x v g) = D2 (f x) (U.map (* f' x) v) hess
  where
    hess = U.zipWith (+) (U.map (* f' x) g) (U.map (* f'' x) (cross v v))

{-# INLINE cross #-}
cross :: U.Vector Double -> U.Vector Double -> U.Vector Double
cross u v = U.concatMap (\dy -> U.map (dy*) u) v

{-# INLINE paramVector2 #-}
paramVector2 :: [Double] -> [AD2]
paramVector2 xs = [ D2 x (U.generate l (\j -> if i == j then 1 else 0)) nil
                  | (i,x) <- zip [0..] xs ]
  where l = length xs ; nil = U.replicate (l*l) 0

