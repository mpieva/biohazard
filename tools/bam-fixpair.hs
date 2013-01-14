{-# LANGUAGE BangPatterns, OverloadedStrings #-}

{- 
How to turn this into a useful validator?

We want both: quickly join separate mates together again, but at the
same time deal with broken files where that doesn't actually work.
Whenever we join mates, we also check if the flags are consistent and
fix them if they aren't.

Quick joining can use a priority queue:  if we have a lone mate and the
mate maps to a later position, we enter it into the pqueue with the
mate's coordinates.  If the mate maps to an earlier position, we can
take it from the pqueue's head.  (Needs some fiddling to work with
multiple reads mapping to the same spot.)

If that fails, it's because a mate is missing or the mate information is
fouled up.  So whenever the head of the pqueue has coordinates in the
past, we take it out and print a warning.  We can still keep this junk
(mem or disk?) or we could discard it.  What hasn't found a mate at the
end can still be discarded or flagged as unpaired.


What to fix:  To make paired records consistent, the MRNM, MPOS and
ISIZE fields should reflect the actual mate's alignment.  This is easy
to fix, but only relevant when that info was wrong to begin with and our
nice pqueue didn't work out.  The flags should also be sanitized,
especially reversed/mate reverses and unmapped/mate unmapped need to go
together, and QC fail should be set for both mates or not at all.

We may need write access to RNAME, POS, MRNM, MPOS, ISIZE, FLAGS on raw
records, too.
-}

import Bio.File.Bam
import Bio.Iteratee
import Bio.PriorityQueue
import Paths_biohazard                          ( version )
import System.IO

import qualified Data.HashMap.Strict as M
import qualified Data.ByteString as S

main :: IO ()
main = addPG (Just version)         >>= \add_pg ->
       enumDefaultInputs >=> run     $
       joinI $ decodeAnyBam          $  \hdr ->
       joinI $ re_pair               $
       joinI $ encodeBamUncompressed (add_pg hdr) $
       mapChunksM_ (liftIO . S.hPut stdout)


-- many more fixes to come or at least checks to come
fixmate :: BamRaw -> BamRaw -> [BamRaw]
fixmate r s | br_isFirstMate r && br_isSecondMate s = [r,s]
            | br_isSecondMate r && br_isFirstMate s = [s,r]
            | otherwise = error "names match, but 1st mate / 2nd mate flags do not"


-- I think this can work with priority queues alone:
--
-- - One contains incomplete pairs ordered by mate position.  When we
--   reach a given position and find the 2nd mate, the minimum in this
--   queue must be the 1st mate (or another 1st mate matching another
--   read we'll find here).
--
-- - One contains incomplete pairs ordered by (hash of) qname.  This one
--   is only used if we missed a mate for some reason.  After we read
--   the whole input, all remaining pairs can be pulled off this queue
--   in order of increasing (hash of) qname.
--
-- - At any given position, we will have a number of 1st mates that have
--   been waiting in the queue and a number of 2nd mates that are coming
--   in from the input.  We dump both sets into a queue by qname, then
--   pull them out in pairs.  Stuff that comes off as anything else than
--   a pair gets queued up again.

re_pair :: MonadIO m => Enumeratee [BamRaw] [BamRaw] m a
{-
re_pair = re_pair' <$> liftIO emptyPQ <*> liftIO emptyPQ
  where
    re_pair' by_pos by_name = eneeCheckIfDone (liftI . none (0::Int))
      where
        -- XXX current state: lone pairs at the end are discarded.  
        none !_ k (EOF mx) = idone (liftI k) (EOF mx)

        none !num !k (Chunk []) = liftI $ none n k

        none !num k (Chunk (r:s:rs)) | br_qname r == br_qname s =
            eneeCheckIfDone (\k' -> none (num+1) k' (Chunk rs)) . k . Chunk $ fixmate r s 

        none !num k (Chunk (r:rs)) = do
            when (num `mod` 0x10000 == 0) $ liftIO $ do s1 <- qsize by_pos
                                                        s2 <- qsize by_name
                                                        hPutStrLn stderr $ "[re-pair] " ++ shows s1 " linked records and "
                                                                           ++ shows s2 " random records queued." 

            top <- getmin by_pos
            case compare (ByPos r) <$> top of
                Nothing ->
                Just LT
                Just EQ ->
                Just GT ->

            case M.lookup (br_qname r) m of
                _ | not (br_isPaired r) -> eneeCheckIfDone (\k' -> go (n+1) m k' (Chunk rs)) . k $ Chunk [r]
                Nothing -> let cp = bamRaw (virt_offset r) (S.copy $ raw_data r)
                           in go (n+1) (M.insert (br_qname cp) cp m) k (Chunk rs)
                Just s  -> eneeCheckIfDone (\k' -> go (n+1) (M.delete (br_qname r) m) k' (Chunk rs)) . k $
                           Chunk $ if br_isFirstMate r then [r,s] else [s,r]
                           -}
