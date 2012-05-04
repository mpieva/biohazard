-- Random useful stuff I didn't know where to put.
module Bio.Util (
    wilson
                ) where

import Data.Number.Erf

-- | calculates the Wilson Score interval.
-- If @(l,m,h) = wilson c x n@, then @m@ is the binary proportion and
-- @(l,h)@ it's @c@-confidence interval for @x@ positive examples out of
-- @n@ observations.  @c@ is typically something like 0.05.

wilson :: Double -> Int -> Int -> (Double, Double, Double)
wilson c x n = ( (m - h) / d, p, (m + h) / d )
  where
    nn = fromIntegral n 
    p  = fromIntegral x / nn

    z = invnormcdf (1-c*0.5)
    h = z * sqrt (( p * (1-p) + 0.25*z*z / nn ) / nn)
    m = p + 0.5 * z * z / nn
    d = 1 + z * z / nn
                
