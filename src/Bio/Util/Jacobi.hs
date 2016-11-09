-- | Jacobi algorithm for Eigen decomposition of symmetric matrices.

module Bio.Util.Jacobi where

import Bio.Prelude

import qualified Data.Vector.Unboxed as U ( Vector, thaw, unsafeFreeze, fromListN, slice )
import qualified Data.Vector.Unboxed.Mutable as U ( read, write, MVector, length, replicate )

type Matrix = U.Vector Double
type Vector = U.Vector Double

type MMatrix s = U.MVector s Double


-- | Decomposes a symmetric matrix into a vector of eigenvalues and a
-- matrix of eigenvectors.

eigen :: Matrix -> ( Vector, [Vector] )
eigen mat = runST (do
    m <- U.thaw mat
    let n = round (sqrt (fromIntegral (U.length m) :: Double))

    v <- U.replicate (n*n) 0
    forM_ [0..n-1] $ \i -> U.write v (n*i+i) 1

    let j_iter !rd = do
            sm <- sum . map abs <$> sequence
                    [ U.read m (p*n+q) | p <- [1..n-2], q <- [p+1..n-1] ]
            case () of
                _ | sm ==  0  -> return ()     -- normal exit: convergence to machine precision
                  | rd == 50  -> return ()     -- 50 iters, weird, but IDGAF.
                  | otherwise -> do
                        let thresh = if rd < 3 then 0.2 * sm / fromIntegral (n*n) else 0
                        forM_ [1..n-1] $ \k ->
                            forM_ [0..k-1] $ \l ->
                                jacobi_rot n thresh (l,k) m v
                        j_iter (rd+1)
    j_iter (0::Int)

    d  <- U.fromListN n <$> forM [0..n-1] (\i -> U.read m (n*i+i))
    v' <- U.unsafeFreeze v
    return (d, [ U.slice i n v' | i <- [0,n..n*n-1] ]))


-- | Performs one Jacobi rotation at @(p,q)@.  We operate on the upper
-- triangle, so at all times @p<q@.  Algorithm inspired by "Numerical
-- Recipes in C: The Art of Scientific Computing" (ISBN 0-521-43108-6)
jacobi_rot :: Int -> Double -> (Int,Int) -> MMatrix s -> MMatrix s -> ST s ()
jacobi_rot n thresh (p,q) m v = do
    a_pp <- U.read m (p*n+p)
    a_pq <- U.read m (p*n+q)
    a_qq <- U.read m (q*n+q)
    let g = 100 * abs a_pq

    case () of
        _ | abs a_pq <= thresh          -> return ()
          | abs a_pp + g == abs a_pp &&
            abs a_qq + g == abs a_qq    -> U.write v (p*n+q) 0
          | otherwise -> do
                let theta = 0.5 * (a_qq-a_pp) / a_pq
                    t     = if theta*theta + 1 == theta*theta
                            then recip $ 2 * theta   -- approx for overflow case
                            else signum theta / ( abs theta + sqrt ( theta*theta +1 ) )
                    c     = recip $ sqrt $ t*t + 1
                    s     = c*t
                    tau   = s / (1+c)

                forM_ [0..n-1] $ \r -> do
                    when (r/=p && r/=q) $ do
                        a_rp <- U.read m $ min (r*n+p) (p*n+r)
                        a_rq <- U.read m $ min (r*n+q) (q*n+r)

                        let a_rp' = a_rp - s * (a_rq + tau * a_rp)
                            a_rq' = a_rq + s * (a_rp - tau * a_rq)

                        U.write m (min (r*n+p) (p*n+r)) a_rp'
                        U.write m (min (r*n+q) (q*n+r)) a_rq'

                    v_rp <- U.read v (p*n+r)
                    v_rq <- U.read v (q*n+r)

                    U.write v (p*n+r) $ v_rp - s * (v_rq + tau * v_rp)
                    U.write v (q*n+r) $ v_rq + s * (v_rp - tau * v_rq)

                U.write m (p*n+q) 0
                U.write m (p*n+p) $ a_pp - t * a_pq
                U.write m (q*n+q) $ a_qq + t * a_pq

