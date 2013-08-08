{-# LANGUAGE BangPatterns, NoMonomorphismRestriction #-}
import Bio.File.Bam
import Bio.Iteratee
import System.Environment
import System.Exit

import qualified Data.ByteString.Char8  as S
import qualified Data.Iteratee          as I

main :: IO ()
main = do
    mq <- getArgs >>= \args -> case (args, reads (head args)) of
            ([ ], _)        -> return 0
            ([_], [(x,[])]) -> return x
            _               -> putStrLn "usage: count-coverage [<min-mapq>]" >> exitFailure

    let putLine nm cv = putStr $ nm ++ '\t' : shows cv "\n"
        printOne refs (r,c) = putLine (S.unpack (sq_name (getRef refs r))) c
        do_count = I.foldl' (\a -> (+) a . snd) 0

    (total,()) <- enumHandle defaultBufSize stdin >=> run $
                  joinI $ decompressBgzf                  $
                  joinI $ decodeBam                       $ \hdr ->
                  joinI $ mapStream decodeBamEntry        $
                  joinI $ filterStream (not . isUnmapped) $
                  joinI $ filterStream ((>= mq) . b_mapq) $
                  joinI $ groupStreamOn b_rname count_cov $        
                  I.zip do_count (I.mapM_ $ printOne $ meta_refs hdr)

    putLine "total" total
        

count_cov :: Monad m => a -> m (Iteratee [BamRec] m Int)
count_cov _ = return $ liftI $ step 0
  where
    step !a (EOF ex) = idone a (EOF ex)
    step !a (Chunk []) = liftI $ step a
    step !a (Chunk (r:rs)) = 
        extend a (b_pos r) (b_pos r + cigarToAlnLen (b_cigar r)) (Chunk rs)

    extend !a !u !v (EOF ex) = idone (a+v-u) (EOF ex)
    extend !a !u !v (Chunk []) = liftI $ extend a u v
    extend !a !u !v (Chunk (r:rs)) 
        | b_pos r <= v = extend a u (max v (b_pos r + cigarToAlnLen (b_cigar r))) (Chunk rs)
        | otherwise    = step (a+v-u) (Chunk (r:rs))




