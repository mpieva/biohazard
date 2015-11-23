{-# LANGUAGE BangPatterns, NoMonomorphismRestriction, FlexibleContexts #-}
import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Bam.Rec
import Bio.Base
import Bio.Iteratee
import Data.Version ( showVersion )
import Paths_biohazard ( version )
import System.Environment
import System.Exit
import System.IO ( hPutStr )

main :: IO ()
main = do
    mq <- getArgs >>= \args -> case (args, reads (head args)) of
            ([ ], _)        -> return (Q 0)
            ([_], [(x,[])]) -> return (Q x)
            _               -> do pn <- getProgName
                                  hPutStr stderr $ pn ++ ", version " ++ showVersion version
                                                ++ "\nUsage: " ++ pn ++ "[<min-mapq>]\n"
                                  exitFailure

    let putLine nm cv = putStr $ nm ++ '\t' : shows cv "\n"

        printOne :: Refs -> (Refseq, Int) -> IO ()
        printOne refs (r,c) = putLine (unpackSeqid (sq_name (getRef refs r))) c

        do_count :: Monad m => Iteratee [(a,Int)] m Int
        do_count = foldStream (\a -> (+) a . snd) 0

    (total,()) <- enumHandle defaultBufSize stdin >=> run                                   $
                  joinI $ decodeAnyBam                                                      $ \hdr ->
                  joinI $ mapMaybeStream ( \br -> case unpackBam br of
                        b | not (isUnmapped b) && b_mapq b >= mq
                            -> Just $! P (b_rname b) (b_pos b) (alignedLength (b_cigar b))
                        _   -> Nothing )                                                    $
                  joinI $ groupStreamOn ref count_cov                                       $
                  zipStreams do_count (mapStreamM_ $ printOne $ meta_refs hdr)

    putLine "total" total

data P = P { ref :: !Refseq, pos :: !Int, alen :: !Int }

count_cov :: Monad m => a -> m (Iteratee [P] m Int)
count_cov _ = return $ liftI $ step 0
  where
    step !a (EOF ex) = idone a (EOF ex)
    step !a (Chunk [    ]) = liftI $ step a
    step !a (Chunk (r:rs)) = extend a (pos r) (pos r + alen r) (Chunk rs)

    extend !a !u !v (EOF ex) = idone (a+v-u) (EOF ex)
    extend !a !u !v (Chunk [    ]) = liftI $ extend a u v
    extend !a !u !v (Chunk (r:rs))
        | pos r <= v = extend a u (max v (pos r + alen r)) (Chunk rs)
        | otherwise  = step (a+v-u) (Chunk (r:rs))






