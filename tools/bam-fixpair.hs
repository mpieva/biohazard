{-# LANGUAGE BangPatterns, OverloadedStrings #-}

{- 
How to turn this into a useful validator?

We want both: quickly join separate mates together again, but at the
same time deal with broken files where that doesn't actually work.
Whenever we join mates, we also check if the flags are consistent and
fix them if they aren't.

In the end, the code will work...

 - splendidly if mates are already adjacent, in which case everything is
   streamed.
 - well if the input is sorted properly, in which case most reads
   stream, but improper pairs need to queue until the mate is reached.
 - reasonably if there are occasional lone mates, which will be queued
   to the very end and sorted by hashed-qname before they are recognized
   and repaired.
 - awkwardly if sorting is violated, flags are wrong or lone mates are
   the rule, because then it degenerates to a full sort by qname.

TODO:
 . actually fix the found pairs
 . actually fix the lone mates
 . deal with consecutive pairs that violate sorting
   (short cut logic:  if consecutive reads form a pair, fix it and pass
   it on; don't fiddle with queues)
 . upgrade to pqueue in external memory
 . better diagnostics
 . useful command line
 . need write access to POS, MRNM, MPOS, ISIZE, FLAGS on raw records
 . a companion that sorts would be cool, but it should be an
   opportunistic sort that is fast on almost sorted files.
 . optional repair:  if 'u' or 'U', but not 'uU' and the mate is
   missing, throw read away  
-}

import Bio.File.Bam
import Bio.Iteratee
import Bio.PriorityQueue
import Data.Binary
import Data.Hashable
import Paths_biohazard                          ( version )
import System.IO
import Text.Printf

import qualified Data.HashMap.Strict as M
import qualified Data.ByteString as S

-- placeholder...
pqconf :: PQ_Conf
pqconf = PQ_Conf 1000 "/var/tmp/"

main :: IO ()
main = addPG (Just version)                               >>= \add_pg ->
       withPQ pqconf                                        $ \right_here ->
       withPQ pqconf                                        $ \in_order   ->
       withPQ pqconf                                        $ \messed_up  ->
       mergeDefaultInputs combineCoordinates >=> run        $ \hdr ->
       re_pair right_here in_order messed_up               =$
       pipeRawBamOutput (add_pg hdr)


-- | Fix a pair of reads.  Right now fixes their order and checks that
-- one is 1st mate, the other 2nd mate.  More fixes to come.
--
-- To do:
-- . set 'unmapped' if missing rname
-- . set 'mate reversed', 'mate unmapped' to value of mate
-- . set 'qc failed' to logical OR of both
-- . unset 'properly paired' if 'unmapped' or 'mate unmapped'
-- . set mrnm to rname of mate
-- . set mpos to pos of mate
-- . set isize

fixmate :: BamRaw -> BamRaw -> [BamRaw]
fixmate r s | br_isFirstMate r && br_isSecondMate s = [r,s]
            | br_isSecondMate r && br_isFirstMate s = [s,r]
            | otherwise = error "names match, but 1st mate / 2nd mate flags do not"

-- | Turns a lone mate into a single.  Right now does nothing, but it
-- should definitely fix some flags.
--
-- To do:
-- . unset 'paired', 'properly paired', '1st mate', '2nd mate', 'mate
--   unmapped', 'mate reversed'
-- . invalidate mrnm, mpos, isize  

divorce :: BamRaw -> BamRaw
divorce b = b

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

-- To ensure proper cleanup, we require the priority ques to be created
-- outside.  Since one is continually reused, it is important that a PQ
-- that is emptied no longer holds on to files on disk.
re_pair :: MonadIO m
        => PQ ByQName -> PQ ByMatePos -> PQ ByQName 
        -> Enumeratee [BamRaw] [BamRaw] m a
re_pair right_here in_order messed_up = eneeCheckIfDone $ go (0::Int) (0::Int)
  where
    note msg = liftIO $ hPutStrLn stderr $ "[re_pair] info:    " ++ msg
    warn msg = liftIO $ hPutStrLn stderr $ "[re_pair] warning: " ++ msg
    err  msg = liftIO $ hPutStrLn stderr $ "[re_pair] error:   " ++ msg

    no_mate_here l b = do warn $ "[" ++ l ++ "] record "
                              ++ shows (br_qname b) (if br_isFirstMate b then "/1" else "/2")
                              ++ " did not have a mate at the right location."
                          let !b' = force_copy b
                          enqueuePQ (byQName b') messed_up

    no_mate_ever b = do err  $ "record " ++ show (br_qname b) ++ " did not have a mate at all."
                        return $ divorce b

    report i o br | (i+o) `mod` 0x20000 /= 0 = return br
                  | otherwise = do hs <- sizePQ right_here
                                   os <- sizePQ in_order
                                   ms <- sizePQ messed_up
                                   let rs = maybe 0 (unRefseq . br_rname) br
                                   let p = maybe 0 br_pos br
                                   note $ printf "@%d/%d, in: %d, out: %d, here: %d, wait: %d, mess: %d" rs p i o hs os ms
                                   return br

    go !num_in !num_out !k = tryHead >>= report num_in num_out >>= go' num_in num_out k
    
    -- At EOF, flush everything.
    go' !num_in !num_out !k Nothing = flush_in_order num_in num_out k

    -- Single read?  Pass through and go on.
    -- Paired read?  Does it belong 'here'?
    go' !num_in !num_out !k (Just r) 
        | not (br_isPaired r) = eneeCheckIfDone (go (num_in+1) (num_out+1)) . k $ Chunk [r]
        | otherwise = peekMinPQ right_here >>= \mm -> case mm of
            -- there's nothing else here, so here becomes redefined
            Nothing             -> enqueue_and_go (num_in+1) num_out k r 

            Just (ByQName _ qq) -> case compare (br_self_pos r) (br_self_pos qq) of
                -- nope, r is out of order and goes to 'messed_up'
                LT -> do warn $ "record " ++ show (br_qname r) ++ " is out of order."
                         let !r' = force_copy r
                         enqueuePQ (byQName r') messed_up 
                         go (num_in+1) num_out k

                -- nope, r comes later.  we need to finish our business here
                GT -> do complete_here (br_self_pos qq)
                         flush_here num_in num_out r k

                -- it belongs here or there is nothing else here
                EQ -> enqueue_and_go (num_in+1) num_out k r 


    -- lonely guy, belongs either here or needs to wait for the mate
    enqueue_and_go num_in num_out k r = do
        if br_self_pos r < br_mate_pos r 
          then let !r' = force_copy r in enqueuePQ (ByMatePos r') in_order
          else enqueuePQ (byQName r) right_here
        go num_in num_out k

    -- Flush the in_order queue to messed_up, since those didn't find
    -- their mate the ordinary way.  Afterwards, flush the messed_up
    -- queue.
    flush_in_order num_in num_out k = getMinPQ in_order >>= \zz -> case zz of
        Just (ByMatePos b) -> do no_mate_here "flush_in_order" b
                                 flush_in_order num_in num_out k
        Nothing -> flush_messed_up num_in num_out k

    -- Flush the messed up queue.  Everything should come off in pairs,
    -- unless something is broken.
    flush_messed_up  num_in num_out k = getMinPQ messed_up >>= flush_mess1 num_in num_out k

    flush_mess1 num_in num_out k Nothing              = return (liftI k)
    flush_mess1 num_in num_out k (Just (ByQName _ a)) = getMinPQ messed_up >>= flush_mess2 num_in num_out k a

    flush_mess2 num_in num_out k a Nothing = do a' <- no_mate_ever a
                                                return (k $ Chunk [a'])
                                                
    flush_mess2 num_in num_out k a b'@(Just (ByQName _ b)) 
        | br_qname a /= br_qname b = do a' <- no_mate_ever a
                                        eneeCheckIfDone (\k' -> flush_mess1 num_in num_out k' b') . k $ Chunk [a']

        | otherwise = eneeCheckIfDone (flush_messed_up num_in (num_out+2)) . k . Chunk $ fixmate a b


    -- Flush the right_here queue.  Everything should come off in pairs,
    -- if not, it goes to messed_up.  When done, loop back to 'go'
    flush_here  num_in num_out r k = getMinPQ right_here >>= flush_here1 num_in num_out r k

    flush_here1 num_in num_out r k Nothing = go' num_in num_out k (Just r)
    flush_here1 num_in num_out r k (Just a) = getMinPQ right_here >>= flush_here2 num_in num_out r k a

    flush_here2 num_in num_out r k (ByQName _ a) Nothing = do no_mate_here "flush_here2/Nothing" a
                                                              flush_here num_in num_out r k
                                                
    flush_here2 num_in num_out r k (ByQName _ a) b'@(Just (ByQName _ b)) 
        | br_qname a /= br_qname b = do no_mate_here "flush_here2/Just" a
                                        flush_here1 num_in num_out r k b'

        | otherwise = eneeCheckIfDone (flush_here num_in (num_out+2) r) . k . Chunk $ fixmate a b


    -- add stuff coming from 'in_order' to 'right_here'
    complete_here pivot = do
            zz <- peekMinPQ in_order
            case zz of
                Nothing -> return ()
                Just (ByMatePos b) 
                       | pivot  > br_mate_pos b -> do dequeuePQ in_order
                                                      no_mate_here "complete_here" b
                                                      complete_here pivot
                    
                       | pivot == br_mate_pos b -> do dequeuePQ in_order
                                                      enqueuePQ (byQName b) right_here
                                                      complete_here pivot

                       | otherwise -> return ()


data ByQName = ByQName !Int !BamRaw

byQName :: BamRaw -> ByQName
byQName b = ByQName (hash $ br_qname b) b

instance Eq ByQName where 
    ByQName ah a == ByQName bh b = 
        (ah, br_qname a) == (bh, br_qname b)

instance Ord ByQName where 
    ByQName ah a `compare` ByQName bh b = 
        (ah, br_qname a) `compare` (bh, br_qname b)

newtype ByMatePos = ByMatePos BamRaw

instance Eq ByMatePos where
    ByMatePos a == ByMatePos b = 
        br_mate_pos a == br_mate_pos b

instance Ord ByMatePos where
    ByMatePos a `compare` ByMatePos b = 
        br_mate_pos a `compare` br_mate_pos b

instance Binary ByQName     -- XXX
instance Binary ByMatePos   -- XXX

instance Sizeable ByQName       -- XXX
instance Sizeable ByMatePos     -- XXX

br_mate_pos :: BamRaw -> (Refseq, Int)
br_mate_pos b = (br_mrnm b, br_mpos b)

br_self_pos :: BamRaw -> (Refseq, Int)
br_self_pos b = (br_rname b, br_pos b)

force_copy :: BamRaw -> BamRaw
force_copy br = bamRaw (virt_offset br) $! S.copy (raw_data br)

