{-# LANGUAGE BangPatterns #-}
import Bio.Bam.Header
import Bio.Bam.Reader
import Bio.Bam.Rec
import Bio.Base
import Bio.Iteratee

main :: IO ()
main = mergeDefaultInputs combineCoordinates >=> run $ \hdr ->
           joinI $ filterStream (not . isUnmapped . unpackBam) $
           joinI $ groupStreamOn (b_rname . unpackBam) (cov_to_wiggle hdr) $
           skipToEof

cov_to_wiggle :: MonadIO m => BamMeta -> Refseq -> m (Iteratee [BamRaw] m ())
cov_to_wiggle hdr rname = return $ liftI step
  where
    step (EOF       mx) = idone () (EOF mx)
    step (Chunk [    ]) = liftI step
    step (Chunk (x:xs)) = do
            let sid = unpackSeqid . sq_name $ meta_refs hdr `getRef` rname
            liftIO $ putStr $ "chrom=" ++ sid ++ " start=" ++ shows (b_pos $ unpackBam x) " step=1\n"
            step' (0::Int) [] (b_pos $ unpackBam x) (Chunk (x:xs))

    step' !cov (e:ends) p           str  | e == p        = step' (cov-1) ends p str

    step' !cov    ends  p (Chunk [    ])                 = liftI (step' cov ends p)
    step' !cov    ends  p (Chunk (x:xs)) | b_pos y == p  = let !e' = b_pos y + alignedLength (b_cigar y)
                                                           in step' (cov+1) (ins e' ends) p (Chunk xs)
        where y = unpackBam x

    step'    _ [      ] _           str                  = step str
    step' !cov    ends  p           str                  = do liftIO $ putStrLn $ show cov
                                                              step' cov ends (p+1) str

    ins a [] = [a]
    ins a (b:bs) | a <= b    = a : b  :  bs
                 | otherwise = b : ins a bs

