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
 . upgrade to pqueue in external memory
 . a companion program that sorts would be cool, but it should be an
   opportunistic sort that is fast on almost sorted files.
-}

import Bio.Base
import Bio.Bam.Header
import Bio.Bam.Raw hiding ( mergeInputs, combineCoordinates )
import Bio.Iteratee
import Bio.PriorityQueue
import Bio.Util                                 ( showNum )
import Control.Monad
import Control.Monad.Trans.Class
import Data.Binary
import Data.Bits
import Data.Hashable
import Data.List
import Paths_biohazard_tools                    ( version )
import System.Console.GetOpt
import System.Environment                       ( getArgs, getProgName )
import System.Exit                              ( exitFailure, exitSuccess )
import System.IO                                ( hPutStrLn )
import Text.Printf

import qualified Data.ByteString as S

data Verbosity = Silent | Errors | Warnings | Notices deriving (Eq, Ord)
data KillMode  = KillNone | KillUu | KillAll deriving (Eq, Ord)

data Config = CF { report_mrnm :: !Bool
                 , report_mpos :: !Bool
                 , report_isize :: !Bool
                 , report_flags :: !Bool
                 , report_fflag :: !Bool
                 , report_ixs :: !Bool
                 , verbosity :: Verbosity
                 , killmode :: KillMode
                 , output :: BamMeta -> Iteratee [BamRaw] IO () }

config0 :: IO Config
config0 = return $ CF True True False True False True Errors KillNone pipeRawBamOutput

options :: [OptDescr (Config -> IO Config)]
options = [
    Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE",
    Option "n" ["dry-run","validate"] (NoArg set_validate) "No output, validate only",
    Option "k" ["kill-lone"]  (NoArg (\c -> return $ c { killmode = KillAll  })) "Delete all lone mates",
    Option "u" ["kill-unmap"] (NoArg (\c -> return $ c { killmode = KillUu   })) "Delete unmapped lone mates",
    Option [ ] ["kill-none"]  (NoArg (\c -> return $ c { killmode = KillNone })) "Never delete lone mates (default)",

    Option "v" ["verbose"]  (NoArg (\c -> return $ c { verbosity = Notices  })) "Print informational messages",
    Option "w" ["warnings"] (NoArg (\c -> return $ c { verbosity = Warnings })) "Print warnings and errors",
    Option [ ] ["errors"]   (NoArg (\c -> return $ c { verbosity = Errors   })) "Print only errors (default)",
    Option "q" ["quiet"]    (NoArg (\c -> return $ c { verbosity = Silent   })) "Print only fatal errors",

    Option "" ["report-mrnm"]  (NoArg (\c -> return $ c { report_mrnm  = True })) "Report wrong mate reference name (default yes)",
    Option "" ["report-mpos"]  (NoArg (\c -> return $ c { report_mpos  = True })) "Report wrong mate position (default yes)",
    Option "" ["report-isize"] (NoArg (\c -> return $ c { report_isize = True })) "Report wrong insert size (default no)",
    Option "" ["report-flags"] (NoArg (\c -> return $ c { report_flags = True })) "Report wrong flags (default yes)",
    Option "" ["report-fflag"] (NoArg (\c -> return $ c { report_fflag = True })) "Report commonly inconsistent flags (default no)",

    Option "" ["no-report-mrnm"]  (NoArg (\c -> return $ c { report_mrnm  = False })) "Do not report wrong mate reference name",
    Option "" ["no-report-mpos"]  (NoArg (\c -> return $ c { report_mpos  = False })) "Do not report wrong mate position",
    Option "" ["no-report-isize"] (NoArg (\c -> return $ c { report_isize = False })) "Do not report wrong insert size",
    Option "" ["no-report-flags"] (NoArg (\c -> return $ c { report_flags = False })) "Do not report wrong flags",
    Option "" ["no-report-fflag"] (NoArg (\c -> return $ c { report_fflag = False })) "Do not report commonly inconsistent flags",
    Option "" ["no-report-fflag"] (NoArg (\c -> return $ c { report_ixs = False })) "Do not report mismatched index fields",

    Option "h?" ["help","usage"] (NoArg usage) "Print this helpful message" ]
  where
    usage _ = do pn <- getProgName
                 let blah = "Usage: " ++ pn ++ " [OPTION...] [FILE...]\n\
                            \Merge BAM files, rearrange them to move mate pairs together, \
                            \output a file with consistent mate pair information."
                 hPutStrLn stderr $ usageInfo blah options
                 exitSuccess

    set_output "-" c = return $ c { output = pipeRawBamOutput }
    set_output  f  c = return $ c { output = writeRawBamFile f }
    set_validate   c = return $ c { output = \_ -> skipToEof }


-- XXX placeholder...
pqconf :: PQ_Conf
pqconf = PQ_Conf 1000 "/var/tmp/"

main :: IO ()
main = do (opts, files, errors) <- getOpt Permute options `fmap` getArgs
          unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
          config <- foldl (>>=) config0 opts
          add_pg <- addPG $ Just version
          withQueues                                           $ \queues ->
            mergeInputs files >=> run                          $ \hdr ->
            re_pair queues config                             =$
            (output config) (add_pg hdr)


-- | Fix a pair of reads.  Right now fixes their order and checks that
-- one is 1st mate, the other 2nd mate.  More fixes to come.

fixmate :: MonadIO m => BamRaw -> BamRaw -> Mating r m [BamRaw]
fixmate r s | br_isFirstMate r && br_isSecondMate s = sequence [go r s, go s r]
            | br_isSecondMate r && br_isFirstMate s = sequence [go s r, go r s]
            | otherwise = liftIO $ do hPutStrLn stderr $ "Names match, but 1st mate / 2nd mate flags do not: " ++ unpackSeqid (br_qname r)
                                      hPutStrLn stderr $ "There is no clear way to fix this file.  Giving up."
                                      exitFailure
  where
    -- position of 5' end
    pos5 a = if br_isReversed a then br_pos a + br_aln_length a else br_pos a

    -- transfer info from b to a
    go a b | null problems = return a
           | otherwise = do infos <- filter (not . null) `fmap` sequence [ m | (_,_,m) <- problems ]
                            unless (null infos) $ liftIO $ hPutStrLn stderr $ message infos
                            return $ mutateBamRaw a $ sequence_ [ m | (_,m,_) <- problems ]
      where
        problems = filter (\(p,_,_) -> not p) checks
        checks = [ (br_mrnm a  == br_rname b,    setMrnm  (br_rname b),  count_mrnm)
                 , (br_mpos a  == br_pos b,      setMpos  (br_pos b),    count_mpos)
                 , (br_isize a == computedIsize, setIsize computedIsize, count_isize)
                 , (br_flag a === computedFlag,  setFlag  computedFlag,  count_flags)
                 , (br_flag a =!= computedFlag,  setFlag  computedFlag,  count_fflag)
                 , (br_indices a == common_indices, setIndices common_indices, count_ixs) ]

        message infos = "fixing " ++ shows (br_qname a `S.append` if br_isFirstMate a then "/1" else "/2")
                        ": \t" ++ intercalate ", " infos
        !computedFlag' = (if br_rname a == invalidRefseq then (.|. flagUnmapped) else id) .
                         (if br_rname b == invalidRefseq then (.|. flagMateUnmapped) else id) .
                         (if br_isReversed b then (.|. flagMateReversed) else (.&. complement flagMateReversed)) .
                         (if br_isUnmapped b then (.|. flagMateUnmapped) else (.&. complement flagMateUnmapped)) .
                         (if br_isFailsQC  b then (.|. flagFailsQC) else id) $
                         br_flag a

        !properly_paired = computedFlag' .&. (flagUnmapped .|. flagMateUnmapped) == 0 && br_rname a == br_rname b
        !computedFlag    = if properly_paired then computedFlag' else computedFlag' .&. complement flagProperlyPaired
        !computedIsize   = if properly_paired then pos5 b - pos5 a else 0

        reduce f | f .&. flagMateUnmapped == 0 = f .&. complement flagFailsQC
                 | otherwise = f .&. complement (flagFailsQC .|. flagMateReversed)

        f1 === f2 = reduce f1 == reduce f2
        f1 =!= f2 = f1 /= f2 && f1 === f2

        onlyIf f m = (\z -> if z then m else "") `fmap` tells f

        count_mrnm  = do modify $ \c -> c { num_mrnm = 1 + num_mrnm c }
                         let ra = unRefseq (br_mrnm a); rb = unRefseq (br_rname b)
                         onlyIf report_mrnm $ printf "MRNM %d is wrong (%d)" ra rb

        count_mpos  = do modify $ \c -> c { num_mpos = 1 + num_mpos c }
                         onlyIf report_mpos $ printf "MPOS %d is wrong (%d)" (br_mpos a) (br_pos b)

        count_isize = do modify $ \c -> c { num_isize = 1 + num_isize c }
                         onlyIf report_isize $ printf "ISIZE %d is wrong (%d)" (br_isize a) computedIsize

        count_flags = do modify $ \c -> c { num_flags = 1 + num_flags c }
                         onlyIf report_flags $ printf "FLAG %03X is wrong (+%03X,-%03X)" (br_flag a) fp fm

        count_fflag = do modify $ \c -> c { num_fflag = 1 + num_fflag c }
                         onlyIf report_fflag $ printf "FLAG %03X is technically wrong (+%03X,-%03X)" (br_flag a) fp fm

        count_ixs = do modify $ \c -> c { num_ixs = 1 + num_ixs c }
                       onlyIf report_ixs $ printf "Index fields %s are wrong (%s)" (show $ br_indices a) (show common_indices)

        fp = computedFlag .&. complement (br_flag a)
        fm = complement computedFlag .&. br_flag a

        index_fields = words "XI XJ YI YJ RG BC"
        br_indices x = [ br_extAsString key x | key <- index_fields ]
        common_indices = zipWith max (br_indices a) (br_indices b)

        setIndices is = do mapM_ removeExt index_fields
                           zipWithM_ appendStringExt index_fields is

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

data MatingStats = MS { total_in   :: !Int
                      , total_out  :: !Int
                      , singletons :: !Int
                      , lone_mates :: !Int
                      , num_mrnm :: !Int
                      , num_mpos :: !Int
                      , num_isize :: !Int
                      , num_flags :: !Int
                      , num_fflag :: !Int
                      , num_ixs :: !Int }

report_stats :: MatingStats -> String
report_stats ms = unlines [
    "number of records read:          " ++ showNum (total_in ms),
    "number of records written:       " ++ showNum (total_out ms),
    "number of true singletons:       " ++ showNum (singletons ms),
    "number of lone mates:            " ++ showNum (lone_mates ms),
    "number of repaired MRNM values:  " ++ showNum (num_mrnm ms),
    "number of repaired MPOS values:  " ++ showNum (num_mpos ms),
    "number of repaired ISIZE values: " ++ showNum (num_isize ms),
    "number of repaired FLAGS values: " ++ showNum (num_flags ms),
    "number of common FLAGS problems: " ++ showNum (num_fflag ms),
    "number of index field problems:  " ++ showNum (num_ixs ms) ]

data Queues = QS { right_here :: !(PQ ByQName)
                 , in_order   :: !(PQ ByMatePos)
                 , messed_up  :: !(PQ ByQName) }

withQueues :: (Queues -> IO r) -> IO r
withQueues k = withPQ pqconf $ \h ->
               withPQ pqconf $ \o   ->
               withPQ pqconf $ \m  ->
               k $ QS h o m

ms0 :: MatingStats
ms0 = MS 0 0 0 0 0 0 0 0 0 0

getSize :: (MonadIO m, Ord a, Binary a, Sizeable a) => (Queues -> PQ a) -> Mating r m Int
getSize sel = getq sel >>= liftIO . sizePQ

enqueue :: (MonadIO m, Ord a, Binary a, Sizeable a) => a -> (Queues -> PQ a) -> Mating r m ()
enqueue a sel = getq sel >>= liftIO . enqueuePQ a

peekMin :: (MonadIO m, Ord a, Binary a, Sizeable a) => (Queues -> PQ a) -> Mating r m (Maybe a)
peekMin sel = getq sel >>= liftIO . peekMinPQ

fetchMin :: (MonadIO m, Ord a, Binary a, Sizeable a) => (Queues -> PQ a) -> Mating r m (Maybe a)
fetchMin sel = getq sel >>= liftIO . getMinPQ

discardMin :: (MonadIO m, Ord a, Binary a, Sizeable a) => (Queues -> PQ a) -> Mating r m ()
discardMin sel = getq sel >>= liftIO . getMinPQ >>= \_ -> return ()


note, warn, err :: MonadIO m => String -> Mating r m ()
note msg = do v <- tells verbosity ; unless (v < Notices)  $ liftIO $ hPutStrLn stderr $ "[fixpair] info:    " ++ msg
warn msg = do v <- tells verbosity ; unless (v < Warnings) $ liftIO $ hPutStrLn stderr $ "[fixpair] warning: " ++ msg
err  msg = do v <- tells verbosity ; unless (v < Errors)   $ liftIO $ hPutStrLn stderr $ "[fixpair] error:   " ++ msg

report' :: MonadIO m => Mating r m ()
report' = do o <- gets total_out
             when (o `mod` 0x40000 == 0) $ do
                     ms <- getSize messed_up
                     note $ printf "out: %d, mess: %d" o ms

report :: MonadIO m => BamRaw -> Mating r m ()
report br = do i <- gets total_in
               o <- gets total_out
               when (i `mod` 0x20000 == 0) $ do
                     hs <- getSize right_here
                     os <- getSize in_order
                     ms <- getSize messed_up
                     let rs = br_rname br ; p = br_pos br
                         at = if rs == invalidRefseq || p == invalidPos
                              then "" else printf "@%d/%d, " (unRefseq rs) p
                     note $ printf "%sin: %d, out: %d, here: %d, wait: %d, mess: %d" (at::String) i o hs os ms

no_mate_here :: MonadIO m => String -> BamRaw -> Mating r m ()
no_mate_here l b = do note $ "[" ++ l ++ "] record "
                          ++ shows (br_qname b) (if br_isFirstMate b then "/1" else "/2")
                          ++ " did not have a mate at the right location."
                      let !b' = force_copy b
                      enqueue (byQName b') messed_up

no_mate_ever :: MonadIO m => BamRaw -> Mating r m ()
no_mate_ever b = do err $ "record " ++ shows (br_qname b) " (" ++
                          shows (br_extAsInt 1 "XI" b) ") did not have a mate at all."
                    modify $ \c -> c { lone_mates = 1 + lone_mates c }
                    kill <- tells killmode
                    case kill of
                        KillAll  -> return ()
                        KillUu   -> unless (br_isUnmapped b) $ yield [divorce b]
                        KillNone -> yield [divorce b]

-- Basically the CPS version of the State Monad.  CPS is necessary to be
-- able to call 'eneeCheckIfDone' in the middle, and that fixes the
-- underlying monad to an 'Iteratee' and the ultimate return type to an
-- 'Iteratee', too.  Pretty to work with, not pretty to look at.
type Sink r m = Stream [BamRaw] -> Iteratee [BamRaw] m r
newtype Mating r m a = Mating { runMating ::
    (a -> MatingStats -> Sink r m -> Queues -> Config -> Iteratee [BamPair] m (Iteratee [BamRaw] m r))
       -> MatingStats -> Sink r m -> Queues -> Config -> Iteratee [BamPair] m (Iteratee [BamRaw] m r) }

instance Monad m => Monad (Mating r m) where
    return a = Mating $ \k -> k a
    m >>=  k = Mating $ \k2 -> runMating m (\a -> runMating (k a) k2)

instance Functor (Mating r m) where
    fmap f m = Mating $ \k -> runMating m (k . f)

instance MonadIO m => MonadIO (Mating r m) where
    liftIO f = Mating $ \k s o q c -> liftIO f >>= \a -> k a s o q c

instance MonadTrans (Mating r) where
    lift m = Mating $ \k s o q c -> lift m >>= \a -> k a s o q c

lift'it :: Monad m => Iteratee [BamPair] m a -> Mating r m a
lift'it m = Mating $ \k s o q c -> m >>= \a -> k a s o q c

tells :: (Config -> a) -> Mating r m a
tells f = Mating $ \k s o q c -> k (f c) s o q c

gets :: (MatingStats -> a) -> Mating r m a
gets f = Mating $ \k s -> k (f s) s

getq :: (Queues -> a) -> Mating r m a
getq f = Mating $ \k s o q -> k (f q) s o q

modify :: (MatingStats -> MatingStats) -> Mating r m ()
modify f = Mating $ \k s -> (k () $! f s)


fetchNext :: MonadIO m => Mating r m (Maybe BamPair)
fetchNext = do r <- lift'it tryHead
               case r of Nothing -> return ()
                         Just (Singleton x) -> do modify $ \s -> s { total_in = 1 + total_in s } ; report x
                         Just (Pair    _ x) -> do modify $ \s -> s { total_in = 2 + total_in s } ; report x
                         Just (LoneMate  x) -> do modify $ \s -> s { total_in = 1 + total_in s } ; report x
               return r

yield :: MonadIO m => [BamRaw] -> Mating r m ()
yield rs = Mating $ \k s o q c -> let !s' = s { total_out = length rs + total_out s }
                                  in eneeCheckIfDone (\o' -> k () s' o' q c) . o $ Chunk rs

-- To ensure proper cleanup, we require the priority queues to be created
-- outside.  Since one is continually reused, it is important that a PQ
-- that is emptied no longer holds on to files on disk.
re_pair :: MonadIO m => Queues -> Config -> Enumeratee [BamPair] [BamRaw] m a
re_pair qs cf = eneeCheckIfDone $ \out -> runMating go finish ms0 out qs cf
   where
    go = fetchNext >>= go'

    -- At EOF, flush everything.
    go' Nothing = peekMin right_here >>= \mm -> case mm of
            Just (ByQName _ _ qq) -> do complete_here (br_self_pos qq)
                                        flush_here Nothing  -- flush_here loops back here
            Nothing               -> flush_in_order  -- this ends the whole operation

    -- Single read?  Pass through and go on.
    -- Paired read?  Does it belong 'here'?
    go' (Just (Singleton x)) = modify (\c -> c { singletons = 1 + singletons c }) >> yield [x] >> go
    go' (Just (Pair    x y)) = fixmate x y >>= yield >> go
    go' (Just (LoneMate  r)) = peekMin right_here >>= \mm -> case mm of

            -- there's nothing else here, so here becomes redefined
            Nothing             -> enqueueThis r >> go

            Just (ByQName _ _ qq) -> case compare (br_self_pos r) (br_self_pos qq) of
                -- nope, r is out of order and goes to 'messed_up'
                LT -> do warn $ "record " ++ show (br_qname r) ++ " is out of order."
                         let !r' = force_copy r
                         enqueue (byQName r') messed_up
                         go

                -- nope, r comes later.  we need to finish our business here
                GT -> do complete_here (br_self_pos qq)
                         flush_here (Just (LoneMate r))

                -- it belongs here or there is nothing else here
                EQ -> enqueueThis r >> go


    -- lonely guy, belongs either here or needs to wait for the mate
    enqueueThis r | br_self_pos r >= br_mate_pos r = enqueue (byQName r) right_here
                  | otherwise             = r' `seq` enqueue (ByMatePos r') in_order
        where r' = force_copy r

    -- Flush the in_order queue to messed_up, since those didn't find
    -- their mate the ordinary way.  Afterwards, flush the messed_up
    -- queue.
    flush_in_order = fetchMin in_order >>= \zz -> case zz of
        Just (ByMatePos b) -> no_mate_here "flush_in_order" b >> flush_in_order
        Nothing            -> flush_messed_up

    -- Flush the messed up queue.  Everything should come off in pairs,
    -- unless something is broken.
    flush_messed_up = fetchMin messed_up >>= flush_mess1

    flush_mess1 Nothing                 = return ()
    flush_mess1 (Just (ByQName _ ai a)) = fetchMin messed_up >>= flush_mess2 ai a

    flush_mess2  _ a Nothing = no_mate_ever a

    flush_mess2 ai a b'@(Just (ByQName _ bi b))
        | ai /= bi || br_qname a /= br_qname b = no_mate_ever a >> report' >> flush_mess1 b'
        | otherwise                            = fixmate a b    >>= yield >> report' >> flush_messed_up


    -- Flush the right_here queue.  Everything should come off in pairs,
    -- if not, it goes to messed_up.  When done, loop back to 'go'
    flush_here  r = fetchMin right_here >>= flush_here1 r

    flush_here1 r Nothing = go' r
    flush_here1 r (Just a) = fetchMin right_here >>= flush_here2 r a

    flush_here2 r (ByQName _ _ a) Nothing = do no_mate_here "flush_here2/Nothing" a
                                               flush_here r

    flush_here2 r (ByQName _ ai a) b'@(Just (ByQName _ bi b))
        | ai /= bi || br_qname a /= br_qname b = no_mate_here "flush_here2/Just" a >> flush_here1 r b'
        | otherwise                            = fixmate a b >>= yield >> flush_here r


    -- add stuff coming from 'in_order' to 'right_here'
    complete_here pivot = do
            zz <- peekMin in_order
            case zz of
                Nothing -> return ()
                Just (ByMatePos b)
                       | pivot  > br_mate_pos b -> do discardMin in_order
                                                      no_mate_here "complete_here" b
                                                      complete_here pivot

                       | pivot == br_mate_pos b -> do discardMin in_order
                                                      enqueue (byQName b) right_here
                                                      complete_here pivot

                       | otherwise -> return ()

    finish () st o _qs _cf = do liftIO $ hPutStrLn stderr $ report_stats st
                                return (liftI o)

data ByQName = ByQName { _bq_hash :: !Int
                       , _bq_alnid :: !Int
                       , _bq_rec :: !BamRaw }

byQName :: BamRaw -> ByQName
byQName b = ByQName (hash $ br_qname b) (br_extAsInt 0 "XI" b) b

instance Eq ByQName where
    ByQName ah ai a == ByQName bh bi b =
        (ah, ai, br_qname a) == (bh, bi, br_qname b)

instance Ord ByQName where
    ByQName ah ai a `compare` ByQName bh bi b =
        (ah, ai, br_qname a) `compare` (bh, bi, br_qname b)

newtype ByMatePos = ByMatePos BamRaw

instance Eq ByMatePos where
    ByMatePos a == ByMatePos b =
        br_mate_pos a == br_mate_pos b

instance Ord ByMatePos where
    ByMatePos a `compare` ByMatePos b =
        br_mate_pos a `compare` br_mate_pos b

instance Binary ByQName where put = undefined ; get = undefined    -- XXX
instance Binary ByMatePos where put = undefined ; get = undefined -- XXX

instance Sizeable ByQName where usedBytes = undefined       -- XXX
instance Sizeable ByMatePos where usedBytes = undefined    -- XXX

br_mate_pos :: BamRaw -> (Refseq, Int)
br_mate_pos b = (br_mrnm b, br_mpos b)

br_self_pos :: BamRaw -> (Refseq, Int)
br_self_pos b = (br_rname b, br_pos b)

force_copy :: BamRaw -> BamRaw
force_copy br = bamRaw (virt_offset br) $! S.copy (raw_data br)



-- | To catch pairs whose mates are adjacent (either because the file
-- has never been sorted or because it has been group-sorted), we apply
-- preprocessing.  The idea is that if we can catch these pairs early,
-- the priority queues never fill up and we save a ton of processing.
-- Now to make the re-pair algorithm work well, we need to merge-sort
-- inputs.  But after that, the pairs have been separated.  So we apply
-- the preprocessing to each input file, then merge then, then run
-- re-pair.

data BamPair = Singleton BamRaw | Pair BamRaw BamRaw | LoneMate BamRaw


mergeInputs :: MonadCatchIO m => [FilePath] -> Enumerator' BamMeta [BamPair] m a
mergeInputs = go0
  where
    go0 [        ] = enumG $ enumHandle defaultBufSize stdin
    go0 (fp0:fps0) = go fp0 fps0

    go fp [       ] = enum1 fp
    go fp (fp1:fps) = mergeEnums' (go fp1 fps) (enum1 fp) combineCoordinates

    enum1 "-" = enumG $ enumHandle defaultBufSize stdin
    enum1  fp = enumG $ enumFile   defaultBufSize    fp

    enumG ee k = ee >=> run $ joinI $ decodeAnyBam $ \h -> quick_pair (k h)


quick_pair :: Monad m => Enumeratee [BamRaw] [BamPair] m a
quick_pair = eneeCheckIfDone go0
  where
    go0 k = tryHead >>= maybe (return $ liftI k) (\x -> go1 x k)

    go1 x k | not (br_isPaired x) = eneeCheckIfDone go0 . k $ Chunk [Singleton x]
            | otherwise           = tryHead >>= maybe (return . k $ Chunk [LoneMate x]) (\y -> go2 x y k)

    go2 x y k | br_qname x == br_qname y = eneeCheckIfDone go0 . k $ Chunk [Pair x y]
              | otherwise                = eneeCheckIfDone (go1 y) . k $ Chunk [LoneMate x]


combineCoordinates :: Monad m => BamMeta -> Enumeratee [BamPair] [BamPair] (Iteratee [BamPair] m) a
combineCoordinates _ = mergeSortStreams (?)
  where u ? v = if (bp_rname u, bp_pos u) < (bp_rname v, bp_pos v) then Less else NotLess

bp_rname :: BamPair -> Refseq
bp_rname (Singleton u) = br_rname u
bp_rname (Pair    u _) = br_rname u
bp_rname (LoneMate  u) = br_rname u

bp_pos :: BamPair -> Int
bp_pos (Singleton u) = br_pos u
bp_pos (Pair    u _) = br_pos u
bp_pos (LoneMate  u) = br_pos u

