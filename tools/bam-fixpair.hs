{-# LANGUAGE BangPatterns, OverloadedStrings #-}

{- 
This is a validator/fixup for paired end BAM files, that is more
efficient than 'samtools sort -n' followed by 'samtools fixmate'.

We want both: to quickly join separate mates together again from the
information about the mate's mapping coordinate, but at the same time
deal with broken files where that doesn't actually work.  Whenever we
join mates, we also check if the flags are consistent and fix them if
they aren't.

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
 . deal with consecutive pairs that violate sorting
   (short cut logic:  if consecutive reads form a pair, fix it and pass
   it on; don't fiddle with queues)
 . upgrade to pqueue in external memory
 . useful command line w/ better control over diagnostics
 . a companion that sorts would be cool, but it should be an
   opportunistic sort that is fast on almost sorted files.
 . optional repair:  if 'u' or 'U', but not 'uU' and the mate is
   missing, throw read away  
-}

import Bio.File.Bam
import Bio.Iteratee
import Bio.PriorityQueue
import Control.Monad
import Data.Binary
import Data.Bits
import Data.Hashable
import Data.List
import Paths_biohazard                          ( version )
import System.IO
import Text.Printf

import qualified Data.HashMap.Strict as M
import qualified Data.ByteString as S

-- XXX placeholder...
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

fixmate :: BamRaw -> BamRaw -> IO [BamRaw]
fixmate r s | br_isFirstMate r && br_isSecondMate s = sequence [go r s, go s r]
            | br_isSecondMate r && br_isFirstMate s = sequence [go s r, go r s]
            | otherwise = error "names match, but 1st mate / 2nd mate flags do not"
  where
    -- position of 5' end
    pos5 a = if br_isReversed a then br_pos a + br_aln_length a else br_pos a

    -- transfer info from b to a
    go a b | null problems = return a
           | otherwise = do unless (null infos) $ liftIO $ hPutStrLn stderr message
                            return $ mutateBamRaw a $ sequence_ [ m | (_,m,_) <- problems ]
      where
        problems = filter (\(p,_,_) -> not p) checks
        checks = [ (br_mrnm a  == br_rname b,    setMrnm  (br_rname b),  printf "MRNM %d is wrong (%d)" ra rb)
                 , (br_mpos a  == br_pos b,      setMpos  (br_pos b),    printf "MPOS %d is wrong (%d)" (br_mpos a) (br_pos b))
                 , (br_isize a == computedIsize, setIsize computedIsize, printf "ISIZE %d is wrong (%d)" (br_isize a) computedIsize)
                 , (br_flag a  == computedFlag,  setFlag  computedFlag,  [] ) ] -- printf "FLAG %03X is wrong (+%03X,-%03X)" (br_flag a) fp fm) ]
          where ra = unRefseq (br_mrnm a); rb = unRefseq (br_rname b)

        message = "fixing " ++ shows (br_qname a `S.append` if br_isFirstMate a then "/1" else "/2") 
                  ": \t" ++ intercalate ", " infos
        infos   = [ m | (_,_,m) <- problems, not (null m) ]
    
        fp = computedFlag .&. complement (br_flag a)
        fm = br_flag a .&. complement computedFlag

        !computedFlag' = (if br_rname a == invalidRefseq then (.|. flagUnmapped) else id) .
                         (if br_rname b == invalidRefseq then (.|. flagMateUnmapped) else id) .
                         (if br_isReversed b then (.|. flagMateReversed) else (.&. complement flagMateReversed)) .
                         (if br_isUnmapped b then (.|. flagMateUnmapped) else (.&. complement flagMateUnmapped)) .
                         (if br_isFailsQC  b then (.|. flagFailsQC) else id) $
                         br_flag a 

        !properly_paired = computedFlag' .&. (flagUnmapped .|. flagMateUnmapped) == 0 && br_rname a == br_rname b
        !computedFlag    = if properly_paired then computedFlag' else computedFlag' .&. complement flagProperlyPaired 
        !computedIsize   = if properly_paired then pos5 b - pos5 a else 0

-- | Turns a lone mate into a single.  Basically removes the pairing
-- related flags and clear the information concerning the mate. 
divorce :: BamRaw -> BamRaw
divorce b = mutateBamRaw b $ do setFlag $ br_flag b .&. complement pair_flags
                                setMrnm invalidRefseq
                                setMpos invalidPos
                                setIsize 0
  where
    pair_flags = flagPaired .|. flagProperlyPaired .|. 
                 flagFirstMate .|. flagSecondMate .|.
                 flagMateUnmapped .|. flagMateReversed

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
                        return [divorce b]

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
    go' !num_in !num_out !k Nothing = peekMinPQ right_here >>= \mm -> case mm of
            Just (ByQName _ qq) -> do complete_here (br_self_pos qq)
                                      flush_here num_in num_out Nothing k -- flush_here loops back here
            Nothing             -> flush_in_order num_in num_out k -- this ends the whole operation

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
                         flush_here num_in num_out (Just r) k

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

    flush_mess2 num_in num_out k a Nothing = no_mate_ever a >>= return . k . Chunk
                                                
    flush_mess2 num_in num_out k a b'@(Just (ByQName _ b)) 
        | br_qname a /= br_qname b = no_mate_ever a >>= eneeCheckIfDone (\k' -> flush_mess1 num_in num_out k' b') . k . Chunk

        | otherwise = liftIO (fixmate a b) >>= eneeCheckIfDone (flush_messed_up num_in (num_out+2)) . k . Chunk


    -- Flush the right_here queue.  Everything should come off in pairs,
    -- if not, it goes to messed_up.  When done, loop back to 'go'
    flush_here  num_in num_out r k = getMinPQ right_here >>= flush_here1 num_in num_out r k

    flush_here1 num_in num_out r k Nothing = go' num_in num_out k r
    flush_here1 num_in num_out r k (Just a) = getMinPQ right_here >>= flush_here2 num_in num_out r k a

    flush_here2 num_in num_out r k (ByQName _ a) Nothing = do no_mate_here "flush_here2/Nothing" a
                                                              flush_here num_in num_out r k
                                                
    flush_here2 num_in num_out r k (ByQName _ a) b'@(Just (ByQName _ b)) 
        | br_qname a /= br_qname b = do no_mate_here "flush_here2/Just" a
                                        flush_here1 num_in num_out r k b'

        | otherwise = liftIO (fixmate a b) >>= eneeCheckIfDone (flush_here num_in (num_out+2) r) . k . Chunk


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

