{-# LANGUAGE BangPatterns #-}
import Bio.File.Bam
import Bio.Iteratee
import System.IO

import qualified Data.Iteratee.ListLike as I

main :: IO ()
main = enumHandle defaultBufSize stdin >=> run $
           joinI $ decodeAnyBam $ \hdr ->
           joinI $ I.mapStream decodeBamEntry $
           joinI $ I.filter (not . isUnmapped) $
           joinI $ groupStreamOn b_rname (cov_to_wiggle hdr) $
           skipToEof

cov_to_wiggle :: MonadIO m => BamMeta -> Refseq -> m (Iteratee [BamRec] m ())
cov_to_wiggle hdr rname = return $ liftI step
  where
    step (EOF mx) = idone () (EOF mx)
    step (Chunk []) = liftI step 
    step (Chunk (x:xs)) = do
            let sid = unpackSeqid . sq_name $ meta_refs hdr `getRef` rname
            liftIO $ putStr $ "chrom=" ++ sid ++ " start=" ++ shows (b_pos x) " step=1\n"
            step' (0::Int) [] (b_pos x) (Chunk (x:xs))

    step' !cov (e:ends) p           str  | e == p       = step' (cov-1) ends p str

    step' !cov    ends  p (Chunk [    ])                = liftI (step' cov ends p)
    step' !cov    ends  p (Chunk (x:xs)) | b_pos x == p = let !e' = b_pos x + cigarToAlnLen (b_cigar x)
                                                          in step' (cov+1) (ins e' ends) p (Chunk xs)

    step'    _ [      ] _           str                 = step str
    step' !cov    ends  p           str                 = do liftIO $ putStrLn $ show cov
                                                             step' cov ends (p+1) str

    ins a [] = [a]
    ins a (a':as) | a <= a'   = (a:a':as)
                  | otherwise = a' : ins a as

