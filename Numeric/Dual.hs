module Numeric.Dual where

-- | "Dual Numbers", or more practically, numbers together with their
-- derivative.  We will have many independent parameters, so the
-- derivative is a vector.

-- 'Dual Double Double' should give a value and a gradient.  'Dual
-- Double (Dual Double Double)' should give value and two derivatives.
-- The Hessian will be encoded as vector of vectors(?)

data Dual a d = D !a !(V.Vector d)

instance (Num a, Num d) => Num (Dual a d) where
    Dual a u + Dual b v = Dual (a+b) (V.zipWith (+) u v)
    Dual a u - Dual b v = Dual (a+b) (V.zipWith (-) u v)
    Dual a u * Dual b v = Dual (a*b) (V.zipWith (+) (V.map (b*) u) (V.map (a*) v))
    negate (Dual a u) = Dual (negate a) (V.map negate v)
    abs d@(Dual a _) = if a >= 0 then d else negate d
    signum (Dual a d) = Dual (signum a) (V.map (const 0) d)
    fromInteger i = Dual i 0 -- oof, not going to work

instance (Fractional a, Fractional d) => Fractional (Dual a d) where
    Dual a u / Dual b v = Dual (a/b) (V.map (w*) $ V.zipWith (-) (V.map (b*) u) (V.map (a*) v))
      where w = recip (b*b)

    recip (Dual a u) = Dual (recip a) (V.map (w*) u)
      where w = recip (a*a)

    fromRational (a :% b) = Dual (fromInteger a / fromInteger b) 0 -- oof, not going to work

