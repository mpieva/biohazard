{-# LANGUAGE BangPatterns, NoMonomorphismRestriction, FlexibleContexts #-}
import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Base
import Bio.Iteratee
import Data.Version ( showVersion )
import Paths_biohazard ( version )
import System.Environment
import System.Exit
import System.IO ( hPutStr )

import qualified Data.Iteratee          as I

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
        printOne refs (r,c) = putLine (unpackSeqid (sq_name (getRef refs r))) c
        do_count = I.foldl' (\a -> (+) a . snd) 0

    (total,()) <- enumHandle defaultBufSize stdin >=> run $
                  joinI $ decodeAnyBam                    $ \hdr ->
                  joinI $ filterStream (\b -> not (br_isUnmapped b) && br_mapq b >= mq) $
                  joinI $ groupStreamOn br_rname count_cov $
                  I.zip do_count (I.mapM_ $ printOne $ meta_refs hdr)

    putLine "total" total


count_cov :: Monad m => a -> m (Iteratee [BamRaw] m Int)
count_cov _ = return $ liftI $ step 0
  where
    step !a (EOF ex) = idone a (EOF ex)
    step !a (Chunk []) = liftI $ step a
    step !a (Chunk (r:rs)) =
        extend a (br_pos r) (br_pos r + br_aln_length r) (Chunk rs)

    extend !a !u !v (EOF ex) = idone (a+v-u) (EOF ex)
    extend !a !u !v (Chunk []) = liftI $ extend a u v
    extend !a !u !v (Chunk (r:rs))
        | br_pos r <= v = extend a u (max v (br_pos r + br_aln_length r)) (Chunk rs)
        | otherwise     = step (a+v-u) (Chunk (r:rs))






