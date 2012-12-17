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
import System.IO

import qualified Data.HashMap.Strict as M
import qualified Data.ByteString as S

main :: IO ()
main = enumDefaultInputs >=> run $
       joinI $ decodeAnyBam $ \hdr ->
       joinI $ re_pair $
       joinI $ encodeBamUncompressed hdr $
       mapChunksM_ (liftIO . S.hPut stdout)


re_pair :: MonadIO m => Enumeratee [BamRaw] [BamRaw] m a
re_pair = eneeCheckIfDone (liftI . go (0::Int) M.empty)
  where
    -- Note: lone pairs left over at the end are discarded.  No harm
    -- done... this is a quick fix anyway.
    go !_ !_ k (EOF mx) = idone (liftI k) (EOF mx)

    go !n !m k (Chunk []) = liftI $ go n m k

    go !n !m k (Chunk (r:s:rs)) | br_qname r == br_qname s =
        eneeCheckIfDone (\k' -> go (n+1) m k' (Chunk rs)) . k . Chunk $
        if br_isFirstMate r then [r,s] else [s,r]

    go !n !m k (Chunk (r:rs)) = do
        if n `mod` 0x10000 == 0 then liftIO $ hPutStrLn stderr $ "[re-pair] " ++ show (M.size m) ++ " records queued."
                              else return ()
        case M.lookup (br_qname r) m of
            _ | not (br_isPaired r) -> eneeCheckIfDone (\k' -> go (n+1) m k' (Chunk rs)) . k $ Chunk [r]
            Nothing -> let cp = bamRaw (virt_offset r) (S.copy $ raw_data r)
                       in go (n+1) (M.insert (br_qname cp) cp m) k (Chunk rs)
            Just s  -> eneeCheckIfDone (\k' -> go (n+1) (M.delete (br_qname r) m) k' (Chunk rs)) . k $
                       Chunk $ if br_isFirstMate r then [r,s] else [s,r]
