{-# LANGUAGE BangPatterns, OverloadedStrings, RecordWildCards, PatternGuards #-}
{-# OPTIONS_GHC -funbox-strict-fields #-}
module Bio.Bam.Index (
    BamIndex,
    readBamIndex,
    readBaiIndex,

    Region(..),
    decodeBamRefseq,
    decodeBamRegions,
    decodeBamUnaligned
) where

import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Iteratee
import Control.Monad
import Data.Array.IArray
import Data.Array.IO                ( IOArray, IOUArray, newArray_, writeArray )
import Data.Array.Unboxed
import Data.Array.Unsafe            ( unsafeFreeze )
import Data.Bits                    ( shiftL, shiftR )
import Data.ByteString              ( ByteString )
import Data.IntMap                  ( IntMap )
import Data.Word                    ( Word64, Word32 )
import System.Directory             ( doesFileExist )
import System.FilePath              ( dropExtension, takeExtension, (<.>) )

import qualified Data.IntMap                    as IM

-- | Full index, unifying BAI and CSI style.  In both cases, we have the
-- binning scheme, parameters are fixed in BAI, but variable in CSI.  BAI
-- has the additional linear index per target sequence.  CSI has
-- `loffset' per bin, which might serve the same purpose(?).

data BamIndex = BamIndex {
    -- | Minshift parameter from CSI
    bi_minshift :: Int,

    -- | Depth parameter from CSI
    bi_depth :: Int,

    -- | Best guess at where the unaligned records start.
    bi_unaln :: !Word64,

    -- | Records for the binning index, where each bin has a list of
    -- chunks belonging to it.
    bi_refseq_bins :: !(Array Refseq Bins),

    -- | Known checkpoints of the form (pos,off) where off is the
    -- virtual offset of the first record crossing pos.
    bi_refseq_ckpoints :: !(Array Refseq Ckpoints) }

-- | Mapping from bin number to vector of clusters.
type Bins = IntMap Chunks
type Chunks = UArray (Int,Edge) Word64

data Edge = Beg | End deriving (Eq, Ord, Ix)



-- | Checkpoints.  Each checkpoint is a position with the virtual offset
-- where the first alignment crossing the position is found.  In BAI, we
-- get this from the 'ioffset' vector, in CSI we get it from the
-- 'loffset' field:  "Given a region [beg,end), we only need to visit
-- chunks whose end file offset is larger than 'ioffset' of the 16kB
-- window containing 'beg'."  (Sounds like a marginal gain, though.)

type Ckpoints = IntMap Word64


-- | Decode only those reads that fall into one of several regions.
-- Strategy:  We will scan the file mostly linearly, but only those
-- regions that are actually needed.  We filter the decoded stuff so
-- that it actually overlaps our regions.
--
-- From the binning index, we get a list of chunks per requested region.
-- Using the checkpoints, we try to eliminate some of the chunks.  The
-- resulting chunks lists are merged, then traversed.  We seek to the
-- beginning of the earliest chunk and start decoding.  Once the virtual
-- file position leaves the chunk, we seek to the next if there is a
-- sufficiently large gap, else we keep going.

data Region = Region Refseq Int Int

decodeBamRegions :: BamIndex -> [Region] ->  Enumeratee Block [BamRaw] m a
decodeBamRegions = undefined
    -- Make a list of chunks from each region, merge and sort them.
    -- We no longer need a chunk when the start coordinate is beyond the
    -- end of the region that generated the chunk.
    -- We have to filter everything for overlapping the region list;
    -- we'll get a colorful mix from the file.


-- | Read any index we can find for a file.  If the file name doesn't
-- have a .bai or .csi extension, we look for the index by adding such
-- an extension and by replacing the extension with these two.  If
-- nothing looks right, we read the file itself anyway.
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
readBaiIndex = i'getString 4 >>= switch
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
            loffset <- endianRead8 LSB
            let key = llim (fromIntegral bin) (3*depth) minshift
            return $! IM.insertWith min key loffset cp

    -- compute left limit of bin
    llim bin dp sf | dp  ==  0 = 0
                   | bin >= ix = (bin - ix) `shiftL` sf
                   | otherwise = llim bin (dp-3) (sf+3)
            where ix = (1 `shiftL` dp - 1) `div` 7

getIndexArrays :: MonadIO m => Int -> Int
               -> (Word32 -> Ckpoints -> Iteratee ByteString m Ckpoints)
               -> ((Ckpoints, Word64) -> Iteratee ByteString m (Ckpoints, Word64))
               -> Iteratee ByteString m BamIndex
getIndexArrays minshift depth addOneCheckpoint addManyCheckpoints = do
    nref <- fromIntegral `liftM` endianRead4 LSB
    if nref < 1
      then return $ BamIndex minshift depth 0
                             (array (toEnum 1, toEnum 0) [])
                             (array (toEnum 1, toEnum 0) [])
      else do
        rbins  <- liftIO ( newArray_ (toEnum 0, toEnum (nref-1)) :: IO (IOArray Refseq Bins) )
        rckpts <- liftIO ( newArray_ (toEnum 0, toEnum (nref-1)) :: IO (IOArray Refseq Ckpoints) )
        mxR <- reduceM 0 nref 0 $ \mx0 r -> do
                nbins <- endianRead4 LSB
                (bins,cpts,mx1) <- reduceM 0 nbins (IM.empty,IM.empty,mx0) $ \(!im,!cp,!mx) _ -> do
                        bin <- endianRead4 LSB -- the "distinct bin"
                        cp' <- addOneCheckpoint bin cp
                        (mx', chunksarr) <- getChunkArray
                        return (IM.insert (fromIntegral bin) chunksarr im, cp', max mx mx')
                (cpts',mx2) <- addManyCheckpoints (cpts,mx1)
                liftIO $ writeArray rbins (toEnum r) bins >> writeArray rckpts (toEnum r) cpts'
                return mx2
        liftM2 (BamIndex minshift depth mxR) (liftIO $ unsafeFreeze rbins) (liftIO $ unsafeFreeze rckpts)

getChunkArray :: MonadIO m => Iteratee ByteString m (Word64, Chunks)
getChunkArray = do
    nchunks <- fromIntegral `liftM` endianRead4 LSB
    chunksarr <- liftIO $ ( newArray_ ((0,Beg), (nchunks-1,End)) :: IO (IOUArray (Int,Edge) Word64) )
    !mx <- reduceM 0 nchunks 0 $ \mx i -> do beg <- endianRead8 LSB
                                             end <- endianRead8 LSB
                                             liftIO $ do writeArray chunksarr (i,Beg) beg
                                                         writeArray chunksarr (i,End) end
                                             return $! max end mx
    (,) mx `liftM` liftIO (unsafeFreeze chunksarr)

{-# INLINE reduceM #-}
reduceM :: (Monad m, Enum ix, Eq ix) => ix -> ix -> a -> (a -> ix -> m a) -> m a
reduceM beg end acc cons = if beg /= end then cons acc beg >>= \n -> reduceM (succ beg) end n cons else return acc


-- | Seek to a given sequence in a Bam file, read those records.  We use
-- the first checkpoint available for the sequence.
-- requires an appropriate index (read separately), and the file must
-- have been opened in such a way as to allow seeking.  Enumerates over
-- the @BamRaw@ records of the correct sequence only, doesn't enumerate
-- at all if the sequence isn't found.

decodeBamRefseq :: Monad m => BamIndex -> Refseq -> Enumeratee [BamRaw] [BamRaw] m a
decodeBamRefseq BamIndex{..} refseq iter
    | bounds bi_refseq_ckpoints `inRange` refseq
    , Just (voff, _) <- IM.minView (bi_refseq_ckpoints ! refseq)
    , voff /= 0 = do seek $ fromIntegral voff
                     breakE ((/=) refseq . br_rname) iter
    | otherwise = return iter

-- | Seek to the part of a Bam file that contains unaligned reads and
-- decode those.  Sort of the dual to @decodeBamRefseq@.  We use the
-- best guess at where the unaligned stuff starts.  If no such guess is
-- available, we decode everything.

decodeBamUnaligned :: Monad m => BamIndex -> Enumeratee [BamRaw] [BamRaw] m a
decodeBamUnaligned BamIndex{..} iter = do when (bi_unaln /= 0) $ seek $ fromIntegral bi_unaln
                                          filterStream (not . isValidRefseq . br_rname) iter

