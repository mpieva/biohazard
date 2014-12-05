{-# LANGUAGE BangPatterns, OverloadedStrings, RecordWildCards, PatternGuards #-}
{-# OPTIONS_GHC -funbox-strict-fields #-}
module Bio.Bam.Index (
    BamIndex,
    readBamIndex,
    readBaiIndex,

    Region(..),
    Subsequence(..),
    eneeBamRefseq,
    eneeBamSubseq,
    eneeBamRegions,
    eneeBamUnaligned
) where

import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Bam.Regions              ( Region(..), Subsequence(..) )
import Bio.Iteratee
import Control.Monad
import Data.Bits                    ( shiftL, shiftR )
import Data.ByteString              ( ByteString )
import Data.Int                     ( Int64 )
import Data.IntMap                  ( IntMap )
import Data.Word                    ( Word32 )
import System.Directory             ( doesFileExist )
import System.FilePath              ( dropExtension, takeExtension, (<.>) )

import qualified Bio.Bam.Regions                as R
import qualified Data.IntMap                    as IM
import qualified Data.Vector                    as V
import qualified Data.Vector.Mutable            as M
import qualified Data.Vector.Unboxed            as U
import qualified Data.Vector.Unboxed.Mutable    as N
import qualified Data.Vector.Algorithms.Intro   as N

-- | Full index, unifying BAI and CSI style.  In both cases, we have the
-- binning scheme, parameters are fixed in BAI, but variable in CSI.
-- Checkpoints are created from the linear index in BAI or from the
-- `loffset' field in CSI.

data BamIndex = BamIndex {
    -- | Minshift parameter from CSI
    minshift :: !Int,

    -- | Depth parameter from CSI
    depth :: !Int,

    -- | Best guess at where the unaligned records start.
    unaln_off :: !Int64,

    -- | Records for the binning index, where each bin has a list of
    -- segments belonging to it.
    refseq_bins :: !(V.Vector Bins),

    -- | Known checkpoints of the form (pos,off) where off is the
    -- virtual offset of the first record crossing pos.
    refseq_ckpoints :: !(V.Vector Ckpoints) }

-- | Mapping from bin number to vector of clusters.
type Bins = IntMap Segments
type Segments = U.Vector (Int64,Int64)


-- | Checkpoints.  Each checkpoint is a position with the virtual offset
-- where the first alignment crossing the position is found.  In BAI, we
-- get this from the 'ioffset' vector, in CSI we get it from the
-- 'loffset' field:  "Given a region [beg,end), we only need to visit
-- chunks whose end file offset is larger than 'ioffset' of the 16kB
-- window containing 'beg'."  (Sounds like a marginal gain, though.)

type Ckpoints = IntMap Int64


-- | Decode only those reads that fall into one of several regions.
-- Strategy:  We will scan the file mostly linearly, but only those
-- regions that are actually needed.  We filter the decoded stuff so
-- that it actually overlaps our regions.
--
-- From the binning index, we get a list of segments per requested
-- region.  Using the checkpoints, we prune them:  if we have a
-- checkpoint to the left of the beginning of the interesting region, we
-- can move the start of each segment forward to the checkpoint.  If
-- that makes the segment empty, it can be droppped.
--
-- The resulting segment lists are merged, then traversed.  We seek to
-- the beginning of the earliest segment and start decoding.  Once the
-- virtual file position leaves the segment or the alignment position
-- moves past the end of the requested region, we move to the next.
-- Moving is a seek if it spans a sufficiently large gap or points
-- backwards, else we just keep going.

-- | A 'Segment' has a start and an end offset, and an "end coordinate"
-- from the originating region.
data Segment = Segment !Int64 !Int64 !Int deriving Show

segmentLists :: BamIndex -> Refseq -> R.Subsequence -> [[Segment]]
segmentLists bi@BamIndex{..} (Refseq ref) (R.Subsequence imap)
        | Just bins <- refseq_bins V.!? fromIntegral ref,
          Just cpts <- refseq_ckpoints V.!? fromIntegral ref
        = [ rgnToSegments bi beg end bins cpts | (beg,end) <- IM.toList imap ]
segmentLists _ _ _ = []

-- from region to list of bins, then to list of segments
rgnToSegments :: BamIndex -> Int -> Int -> Bins -> Ckpoints -> [Segment]
rgnToSegments bi@BamIndex{..} beg end bins cpts =
    [ Segment boff' eoff end
    | bin <- binList bi beg end
    , (boff,eoff) <- maybe [] U.toList $ IM.lookup bin bins
    , let boff' = max boff cpt
    , boff' < eoff ]
  where
    !cpt = maybe 0 snd $ lookupLE beg cpts

-- list of bins for given range of coordinates, from Heng's horrible code
binList :: BamIndex -> Int -> Int -> [Int]
binList BamIndex{..} beg end = binlist' 0 (minshift + 3*depth) 0
  where
    binlist' l s t = if l > depth then [] else [b..e] ++ loop
      where
        b = t + beg `shiftR` s
        e = t + (end-1) `shiftR` s
        loop = binlist' (l+1) (s-3) (t + 1 `shiftL` (3*l))


-- | Merges two lists of segments.  Lists must be sorted, the merge sort
-- merges overlapping segments into one.
infix 4 ~~
(~~) :: [Segment] -> [Segment] -> [Segment]
Segment a b e : xs ~~ Segment u v f : ys
    |          b < u = Segment a b e : (xs ~~ Segment u v f : ys)     -- no overlap
    | a < u && b < v = Segment a v (max e f) : (xs ~~ ys)             -- some overlap
    |          b < v = Segment u v (max e f) : (xs ~~ ys)             -- contained
    | v < a          = Segment u v f : (xs ~~ Segment a b e : ys)     -- no overlap
    | u < a          = Segment u b (max e f) : (xs ~~ ys)             -- some overlap
    | otherwise      = Segment a b (max e f) : (xs ~~ ys)             -- contained
[] ~~ ys = ys
xs ~~ [] = xs


-- | Reads any index we can find for a file.  If the file name has a
-- .bai or .csi extension, we read it.  Else we look for the index by
-- adding such an extension and by replacing the extension with these
-- two, and finally in the file itself.  The first file that exists and
-- can actually be parsed, is used.
readBamIndex :: FilePath -> IO BamIndex
readBamIndex fp | takeExtension fp == ".bai" = fileDriver readBaiIndex fp
                | takeExtension fp == ".csi" = fileDriver readBaiIndex fp
                | otherwise = try               (fp <.> "bai") $
                              try (dropExtension fp <.> "bai") $
                              try               (fp <.> "csi") $
                              try (dropExtension fp <.> "csi") $
                              fileDriver readBaiIndex fp
  where
    try f k = do e <- doesFileExist f
                 if e then do r <- enumFile defaultBufSize f readBaiIndex >>= tryRun
                              case r of Right                     ix -> return ix
                                        Left (IterStringException _) -> k
                      else k

-- Read in dex in BAI or CSI format, recognized automatically.
readBaiIndex :: MonadIO m => Iteratee ByteString m BamIndex
readBaiIndex = iGetString 4 >>= switch
  where
    switch "BAI\1" = getIndexArrays 14 5 (const return) getIntervals
    switch "CSI\1" = do minshift <- fromIntegral `liftM` endianRead4 LSB
                        depth <- fromIntegral `liftM` endianRead4 LSB
                        endianRead4 LSB >>= dropStream . fromIntegral -- aux data
                        getIndexArrays minshift depth (addOneCheckpoint minshift depth) return
    switch magic   = throwErr . iterStrExc $ "index signature " ++ show magic ++ " not recognized"

    -- Read the intervals.  Each one becomes a checkpoint.
    getIntervals (cp,mx0) = do
        nintv <- fromIntegral `liftM` endianRead4 LSB
        reduceM 0 nintv (cp,mx0) $ \(!im,!mx) int -> do
            oo <- fromIntegral `liftM` endianRead8 LSB
            return (if oo == 0 then im else IM.insert (int * 0x4000) oo im, max mx oo)

    -- Insert one checkpoint.  If we already have an entry (can happen
    -- if it comes from a different bin), we conservatively take the min
    addOneCheckpoint minshift depth bin cp = do
            loffset <- fromIntegral `liftM` endianRead8 LSB
            let key = llim (fromIntegral bin) (3*depth) minshift
            return $! IM.insertWith min key loffset cp

    -- compute left limit of bin
    llim bin dp sf | dp  ==  0 = 0
                   | bin >= ix = (bin - ix) `shiftL` sf
                   | otherwise = llim bin (dp-3) (sf+3)
            where ix = (1 `shiftL` dp - 1) `div` 7

getIndexArrays :: MonadIO m => Int -> Int
               -> (Word32 -> Ckpoints -> Iteratee ByteString m Ckpoints)
               -> ((Ckpoints, Int64) -> Iteratee ByteString m (Ckpoints, Int64))
               -> Iteratee ByteString m BamIndex
getIndexArrays minshift depth addOneCheckpoint addManyCheckpoints = do
    nref <- fromIntegral `liftM` endianRead4 LSB
    if nref < 1
      then return $ BamIndex minshift depth 0 V.empty V.empty
      else do
        rbins  <- liftIO $ M.new nref
        rckpts <- liftIO $ M.new nref
        mxR <- reduceM 0 nref 0 $ \mx0 r -> do
                nbins <- endianRead4 LSB
                (bins,cpts,mx1) <- reduceM 0 nbins (IM.empty,IM.empty,mx0) $ \(!im,!cp,!mx) _ -> do
                        bin <- endianRead4 LSB -- the "distinct bin"
                        cp' <- addOneCheckpoint bin cp
                        segsarr <- getSegmentArray
                        let !mx' = if U.null segsarr then mx else max mx (snd (U.last segsarr))
                        return (IM.insert (fromIntegral bin) segsarr im, cp', mx')
                (cpts',mx2) <- addManyCheckpoints (cpts,mx1)
                liftIO $ M.write rbins r bins >> M.write rckpts r cpts'
                return mx2
        liftM2 (BamIndex minshift depth mxR) (liftIO $ V.unsafeFreeze rbins) (liftIO $ V.unsafeFreeze rckpts)

-- | Reads the list of segments from an index file and makes sure
-- it is sorted.
getSegmentArray :: MonadIO m => Iteratee ByteString m Segments
getSegmentArray = do
    nsegs <- fromIntegral `liftM` endianRead4 LSB
    segsarr <- liftIO $ N.new nsegs
    loopM 0 nsegs $ \i -> do beg <- fromIntegral `liftM` endianRead8 LSB
                             end <- fromIntegral `liftM` endianRead8 LSB
                             liftIO $ N.write segsarr i (beg,end)
    liftIO $ N.sort segsarr >> U.unsafeFreeze segsarr

{-# INLINE reduceM #-}
reduceM :: (Monad m, Enum ix, Eq ix) => ix -> ix -> a -> (a -> ix -> m a) -> m a
reduceM beg end acc cons = if beg /= end then cons acc beg >>= \n -> reduceM (succ beg) end n cons else return acc

{-# INLINE loopM #-}
loopM :: (Monad m, Enum ix, Eq ix) => ix -> ix -> (ix -> m ()) -> m ()
loopM beg end k = if beg /= end then k beg >> loopM (succ beg) end k else return ()


-- | Seeks to a given sequence in a Bam file and enumerates only those
-- records aligning to that reference.  We use the first checkpoint
-- available for the sequence.  This requires an appropriate index, and
-- the file must have been opened in such a way as to allow seeking.
-- Enumerates over the @BamRaw@ records of the correct sequence only,
-- doesn't enumerate at all if the sequence isn't found.

eneeBamRefseq :: Monad m => BamIndex -> Refseq -> Enumeratee [BamRaw] [BamRaw] m a
eneeBamRefseq BamIndex{..} (Refseq r) iter
    | Just ckpts <- refseq_ckpoints V.!? fromIntegral r
    , Just (voff, _) <- IM.minView ckpts
    , voff /= 0 = do seek $ fromIntegral voff
                     breakE ((Refseq r /=) . br_rname) iter
    | otherwise = return iter

-- | Seeks to the part of a Bam file that contains unaligned reads and
-- enumerates those.  Sort of the dual to 'eneeBamRefseq'.  We use the
-- best guess at where the unaligned stuff starts.  If no such guess is
-- available, we decode everything.

eneeBamUnaligned :: Monad m => BamIndex -> Enumeratee [BamRaw] [BamRaw] m a
eneeBamUnaligned BamIndex{..} iter = do when (unaln_off /= 0) $ seek $ fromIntegral unaln_off
                                        filterStream (not . isValidRefseq . br_rname) iter

-- | Enumerates one 'Segment'.  Seeks to the start offset, unless
-- reading over the skipped part looks cheaper.  Enumerates until we
-- either cross the end offset or the max position.
eneeBamSegment :: Monad m => Segment -> Enumeratee [BamRaw] [BamRaw] m r
eneeBamSegment (Segment beg end mpos) out = do
    -- seek if it's a backwards seek or more than 512k forwards
    peekStream >>= \x -> case x of
        Just br | beg <= o && beg + 0x8000 > o -> return ()
            where o = fromIntegral $ virt_offset br
        _                                      -> seek $ fromIntegral beg

    let in_segment br = virt_offset br <= fromIntegral end && br_pos br <= mpos
    takeWhileE in_segment out

eneeBamSubseq :: Monad m => BamIndex -> Refseq -> R.Subsequence -> Enumeratee [BamRaw] [BamRaw] m a
eneeBamSubseq bi ref subs
    = let segs = foldr (~~) [] $ segmentLists bi ref subs
          olap br = br_rname br == ref && R.overlaps (br_pos br) (br_pos br + br_aln_length br) subs
      in foldr ((>=>) . eneeBamSegment) return segs ><> filterStream olap

eneeBamRegions :: Monad m => BamIndex -> [R.Region] -> Enumeratee [BamRaw] [BamRaw] m a
eneeBamRegions bi = foldr ((>=>) . uncurry (eneeBamSubseq bi)) return . R.toList . R.fromList


lookupLE :: IM.Key -> IM.IntMap a -> Maybe (IM.Key, a)
lookupLE k m = case ma of
	Just a               -> Just (k,a)
	Nothing | IM.null m1 -> Nothing
                | otherwise  -> Just $ IM.findMax m1
  where (m1,ma,_) = IM.splitLookup k m
