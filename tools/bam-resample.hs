{-# LANGUAGE BangPatterns #-}
-- Resample m out of n `virtual' BAM records.
--
-- Strategy for fair down sampling:  we first count the number of
-- records, then scan again to sample.  Input must be grouped by QNAME
-- (sorted by QNAME is fine).
--
-- Usage: resample [NUM] [FILE...]

import Bio.File.Bam
import Bio.Iteratee
import Paths_biohazard ( version )
import System.Environment
import System.IO
import System.Random

import qualified Data.ByteString as S
import qualified Data.Iteratee as I

main :: IO ()
main = do
    num_ : files <- getArgs
    num <- readIO num_

    hPutStr stderr "counting... "
    total <- enumInputs files >=> run $
             joinI $ decodeAnyBam $ \_hdr ->
             joinI $ groupOn br_qname $
             foldStream (\a _ -> 1+a) 0
    hPutStr stderr $ shows total " records.\n"

    add_pg <- addPG (Just version)
    enumInputs files >=> run $
             joinI $ decodeAnyBam $ \hdr ->
             joinI $ groupOn br_qname $
             joinI $ resample num total $
             joinI $ encodeBam (add_pg hdr) $
             mapChunksM_ (S.hPut stdout)


resample :: MonadIO m => Int -> Int -> Enumeratee [[BamRaw]] [BamRaw] m a
resample m0 n0 | m0 > n0 = error "upsampling requested"
resample m0 n0 = eneeCheckIfDone (go m0 n0)
  where
    go  !m !n k = I.tryHead >>= maybe (return (liftI k)) (go' m n k)
    go' !m !n k a = do r <- liftIO $ randomRIO (0,n-1)
                       if r < m
                         then eneeCheckIfDone (go (m-1) (n-1)) . k $ Chunk a
                         else go m (n-1) k
    
groupOn :: (Monad m, Eq b) => (a -> b) -> Enumeratee [a] [[a]] m c
groupOn f = eneeCheckIfDone (\k -> I.tryHead >>= maybe (return $ liftI k) (\a -> go k [a] (f a)))
  where
    go  k acc fa = I.tryHead >>= maybe (return . k $ Chunk [reverse acc]) (go' k acc fa)
    go' k acc fa b | fa == f b = go k (b:acc) fa 
                   | otherwise = eneeCheckIfDone (\k' -> go k' [b] (f b)) . k $ Chunk [reverse acc]