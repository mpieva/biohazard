module Bio.Util (
    wilson, invnormcdf, choose,
    estimateComplexity, showNum, showOOM,
    log1p, expm1, (<#>),
    lsum, llerp,
    sigmoid2, isigmoid2
                ) where

import Data.List ( foldl1' )
import Data.Char ( intToDigit )

-- ^ Random useful stuff I didn't know where to put.

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

showNum :: Show a => a -> String
showNum = triplets [] . reverse . show
  where
    triplets acc [] = acc
    triplets acc (a:[]) = a:acc
    triplets acc (a:b:[]) = b:a:acc
    triplets acc (a:b:c:[]) = c:b:a:acc
    triplets acc (a:b:c:s) = triplets (',':c:b:a:acc) s

showOOM :: Double -> String
showOOM x | x < 0 = '-' : showOOM (negate x)
          | otherwise = findSuffix (x*10) ".kMGTPEZY"
  where
    findSuffix _ [] = "many"
    findSuffix y (s:ss) | y < 100  = intToDigit (round y `div` 10) : case (round y `mod` 10, s) of
                                            (0,'.') -> [] ; (0,_) -> [s] ; (d,_) -> [s, intToDigit d]
                        | y < 1000 = intToDigit (round y `div` 100) : intToDigit ((round y `mod` 100) `div` 10) :
                                            if s == '.' then [] else [s]
                        | y < 10000 = intToDigit (round y `div` 1000) : intToDigit ((round y `mod` 1000) `div` 100) :
                                            '0' : if s == '.' then [] else [s]
                        | otherwise = findSuffix (y*0.001) ss

-- Stolen from Lennart Augustsson's erf package, who in turn took it rom
-- http://home.online.no/~pjacklam/notes/invnorm/ Accurate to about 1e-9.
invnormcdf :: (Ord a, Floating a) => a -> a
invnormcdf p =
    let a1 = -3.969683028665376e+01
        a2 =  2.209460984245205e+02
        a3 = -2.759285104469687e+02
        a4 =  1.383577518672690e+02
        a5 = -3.066479806614716e+01
        a6 =  2.506628277459239e+00

        b1 = -5.447609879822406e+01
        b2 =  1.615858368580409e+02
        b3 = -1.556989798598866e+02
        b4 =  6.680131188771972e+01
        b5 = -1.328068155288572e+01

        c1 = -7.784894002430293e-03
        c2 = -3.223964580411365e-01
        c3 = -2.400758277161838e+00
        c4 = -2.549732539343734e+00
        c5 =  4.374664141464968e+00
        c6 =  2.938163982698783e+00

        d1 =  7.784695709041462e-03
        d2 =  3.224671290700398e-01
        d3 =  2.445134137142996e+00
        d4 =  3.754408661907416e+00

        pLow = 0.02425

        nan = 0/0

    in  if p < 0 then
            nan
        else if p == 0 then
            -1/0
        else if p < pLow then
            let q = sqrt(-2 * log p)
            in  (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
                 ((((d1*q+d2)*q+d3)*q+d4)*q+1)
        else if p < 1 - pLow then
            let q = p - 0.5
                r = q*q
            in  (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q /
                (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
        else if p <= 1 then
            - invnormcdf (1 - p)
        else
            nan


-- | Try to estimate complexity of a whole from a sample.  Suppose we
-- sampled @total@ things and among those @singles@ occured only once.
-- How many different things are there?
--
-- Let the total number be @m@.  The copy number follows a Poisson
-- distribution with paramter @\lambda@.  Let @z := e^{\lambda}@, then
-- we have:
--
--   P( 0 ) = e^{-\lambda} = 1/z
--   P( 1 ) = \lambda e^{-\lambda} = ln z / z
--   P(>=1) = 1 - e^{-\lambda} = 1 - 1/z
--
--   singles = m ln z / z
--   total   = m (1 - 1/z)
--
--   D := total/singles = (1 - 1/z) * z / ln z
--   f := z - 1 - D ln z = 0
--
-- To get @z@, we solve using Newton iteration and then substitute to
-- get @m@:
--
--   df/dz = 1 - D/z
--   z' := z - z (z - 1 - D ln z) / (z - D)
--   m = singles * z /log z
--
-- It converges as long as the initial @z@ is large enough, and @10D@
-- (in the line for @zz@ below) appears to work well.

estimateComplexity :: (Integral a, Floating b, Ord b) => a -> a -> Maybe b
estimateComplexity total singles | total   <= singles = Nothing
                                 | singles <= 0       = Nothing
                                 | otherwise          = Just m
  where
    d = fromIntegral total / fromIntegral singles
    step z = z * (z - 1 - d * log z) / (z - d)
    iter z = case step z of zd | abs zd < 1e-12 -> z
                               | otherwise -> iter $! z-zd
    zz = iter $! 10*d
    m = fromIntegral singles * zz / log zz


-- | Computes @log (exp x + exp y)@ without leaving the log domain and
-- hence without losing precision.
infixl 5 <#>
{-# INLINE (<#>) #-}
(<#>) :: (Floating a, Ord a) => a -> a -> a
x <#> y = if x >= y then x + log1p (exp (y-x)) else y + log1p (exp (x-y))

-- | Computes @log (1+x)@ to a relative precision of @10^-8@ even for
-- very small @x@.  Stolen from http://www.johndcook.com/cpp_log_one_plus_x.html
{-# INLINE log1p #-}
log1p :: (Floating a, Ord a) => a -> a
log1p x | x < -1 = error "log1p: argument must be greater than -1"
        -- x is large enough that the obvious evaluation is OK:
        | x > 0.0001 || x < -0.0001 = log $ 1 + x
        -- Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
        -- Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8:
        | otherwise = (1 - 0.5*x) * x


-- | Computes @exp x - 1@ to a relative precision of @10^-10@ even for
-- very small @x@.  Stolen from http://www.johndcook.com/cpp_expm1.html
expm1 :: (Floating a, Ord a) => a -> a
expm1 x | x > -0.00001 && x < 0.00001 = (1 + 0.5 * x) * x       -- Taylor approx
        | otherwise                   = exp x - 1               -- direct eval

-- | Computes \( \log ( \sum_i e^{x_i} ) \) sensibly.  The list must be
-- sorted in descending(!) order.
{-# INLINE lsum #-}
lsum :: (Floating a, Ord a) => [a] -> a
lsum xs = foldl1' (\x y -> if x >= y then x + log1p (exp (y-x)) else err) xs
    where err = error $ "lsum: argument list must be in descending order"

-- | Computes \( \log \left( c e^x + (1-c) e^y \right) \).
{-# INLINE llerp #-}
llerp :: (Floating a, Ord a) => a -> a -> a -> a
llerp c x y | c <= 0.0  = y
            | c >= 1.0  = x
            | x >= y    = log     c  + x + log1p ( (1-c)/c * exp (y-x) )
            | otherwise = log1p (-c) + y + log1p ( c/(1-c) * exp (x-y) )

-- | Binomial coefficient:  @n `choose` k == n! / ((n-k)! k!)@
{-# INLINE choose #-}
choose :: Integral a => a -> a -> a
n `choose` k = product [n-k+1 .. n] `div` product [2..k]


-- | Kind-of sigmoid function that maps the reals to the interval
-- @[0,1)@.  Good to compute a probability without introducing boundary
-- conditions.
sigmoid2 :: (Num a, Fractional a, Floating a) => a -> a
sigmoid2 x = y*y where y = (exp x - 1) / (exp x + 1)

-- | Inverse of 'sigmoid2'.
isigmoid2 :: (Num a, Fractional a, Floating a) => a -> a
isigmoid2 y = log $ (1 + sqrt y) / (1 - sqrt y)


