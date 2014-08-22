module Bio.Genocall.Matrix where

import Bio.Base
import qualified Data.Vector.Unboxed    as V

-- | 1x4 vector used to represent mixes of nucleotide (e.g. genotypes)
newtype Vec = Vec { unVec :: V.Vector Double } deriving Show

-- | 4x4 matrix used for substitution models
newtype Matrix = Matrix { unMatrix :: V.Vector Double }

data To = Nucleotide :-> Nucleotide

infix 9 :->
infix 8 !

(!) :: Matrix -> To -> Double
(!) (Matrix m) (N x :-> N y) = m V.! fromIntegral (4*y+x)

construct :: (To -> Double) -> Matrix
construct f = Matrix $ V.fromListN 16 [ f $ N x :-> N y | y <- [0..3], x <- [0..3] ]

fromList :: [Double] -> Matrix
fromList = Matrix . V.fromListN 16

zeroMatrix :: Matrix
zeroMatrix = Matrix $ V.replicate 16 0

dot :: Vec -> Vec -> Double
dot (Vec x) (Vec y) = V.sum $ V.zipWith (*) x y

mult :: Matrix -> Vec -> Vec
mult (Matrix m) (Vec v) = Vec . V.generate 4 $ \y -> V.sum $ V.zipWith (*) (V.slice (4*y) 4 m) v

update :: Matrix -> [(To, Double)] -> Matrix
update (Matrix m) ps = Matrix $ m V.// [(fromIntegral (4*y+x), a) | (N x :-> N y, a) <- ps ]
