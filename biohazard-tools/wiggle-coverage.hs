{-# LANGUAGE BangPatterns #-}
import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Base
import Bio.Iteratee

main :: IO ()
main = mergeDefaultInputs combineCoordinates >=> run $ \hdr ->
           joinI $ filterStream (not . br_isUnmapped) $
           joinI $ groupStreamOn br_rname (cov_to_wiggle hdr) $
           skipToEof

cov_to_wiggle :: MonadIO m => BamMeta -> Refseq -> m (Iteratee [BamRaw] m ())
cov_to_wiggle hdr rname = return $ liftI step
  where
    step (EOF       mx) = idone () (EOF mx)
    step (Chunk [    ]) = liftI step
    step (Chunk (x:xs)) = do
            let sid = unpackSeqid . sq_name $ meta_refs hdr `getRef` rname
            liftIO $ putStr $ "chrom=" ++ sid ++ " start=" ++ shows (br_pos x) " step=1\n"
            step' (0::Int) [] (br_pos x) (Chunk (x:xs))

    step' !cov (e:ends) p           str  | e == p        = step' (cov-1) ends p str

    step' !cov    ends  p (Chunk [    ])                 = liftI (step' cov ends p)
    step' !cov    ends  p (Chunk (x:xs)) | br_pos x == p = let !e' = br_pos x + br_aln_length x
                                                           in step' (cov+1) (ins e' ends) p (Chunk xs)

    step'    _ [      ] _           str                  = step str
    step' !cov    ends  p           str                  = do liftIO $ putStrLn $ show cov
                                                              step' cov ends (p+1) str

    ins a [] = [a]
    ins a (b:bs) | a <= b    = a : b  :  bs
                 | otherwise = b : ins a bs

