-- | Quality filters adapted from old pipeline.
--
-- TODO: - "SAGE" filter (enforce 17nt reads)?
--       - "rim job" (try to detect glass border)?

{-# LANGUAGE BangPatterns #-}
module Bio.File.Bam.Filter (
    filterPairs, QualFilter,
    complexSimple, complexEntropy,
    qualityAverage, qualityMinimum,
    qualityFromOldIllumina, qualityFromNewIllumina
                           ) where

import Bio.File.Bam
import Bio.Iteratee
import Data.Bits
import Data.Word ( Word8 )

import qualified Data.ByteString     as S
import qualified Data.Iteratee       as I
import qualified Data.Vector.Generic as V

-- | A filter/transformation applied to pairs of reads.  We supply a
-- predicate to be applied to single reads and one to be applied to
-- pairs, tha latter can get incomplete pairs, too, if mates have been
-- separated or filtered asymmetrically.

filterPairs :: Monad m => (BamRec -> [BamRec])
                       -> (Maybe BamRec -> Maybe BamRec -> [BamRec])
                       -> Enumeratee [BamRec] [BamRec] m a
filterPairs ps pp = eneeCheckIfDone step
  where
    step k = I.tryHead >>= step' k
    step' k Nothing = return $ liftI k
    step' k (Just b)
        | isPaired b = I.tryHead >>= step'' k b
        | otherwise  = case ps b of [] -> step k ; b' -> eneeCheckIfDone step . k $ Chunk b'

    step'' k b Nothing = case pp (Just b) Nothing of 
                            [] -> return $ liftI k
                            b' -> return $ k $ Chunk b'

    step'' k b (Just c)
        | b_rname b /= b_rname c || not (isPaired c) =
                let b' = if isSecondMate b then pp Nothing (Just b) else pp (Just b) Nothing
                in case b' of [] -> step' k (Just c)
                              _  -> eneeCheckIfDone (\k' -> step' k' (Just c)) . k $ Chunk b'

        | isFirstMate c && isSecondMate b = step''' k c b
        | otherwise                       = step''' k b c

    step''' k b c = case pp (Just b) (Just c) of [] -> step k
                                                 b' -> eneeCheckIfDone step . k $ Chunk b'


-- | A quality filter is simply a transformation on @BamRec@s.  By
-- convention, quality filters should set @flagFailsQC@, a further step
-- can then remove the failed reads.  Filtering of individual reads
-- tends to result in mate pairs with inconsistent flags, which in turn
-- will result in lone mates and all sort of troubles with programs that
-- expect non-broken BAM files.  It is therefore recommended to use
-- @pairFilter@ with suitable predicates to do the post processing.

type QualFilter = BamRec -> BamRec

{-# INLINE count #-}
count :: (V.Vector v a, Eq a) => a -> v a -> Int
count x v = V.foldl' (\acc y -> if x == y then acc+1 else acc) 0 v

-- | Simple complexity filter aka "Nancy Filter".  A read is considered
-- not-sufficiently-complex if the most common base accounts for greater
-- than the @cutoff@ fraction of all non-N bases.
{-# INLINE complexSimple #-}
complexSimple :: Double -> QualFilter
complexSimple r b = if p then b else b'
  where
    b' = setQualFlag 'C' $ b { b_flag = b_flag b .|. flagFailsQC }
    p  = let counts = [ count x $ b_seq b | x <- properBases ]
             lim = floor $ r * fromIntegral (sum counts)
         in all (<= lim) counts

-- | Filter on order zero empirical entropy.  Entropy per base must be
-- greater than cutoff.
{-# INLINE complexEntropy #-}
complexEntropy :: Double -> QualFilter
complexEntropy r b = if p then b else b'
  where
    b' = setQualFlag 'C' $ b { b_flag = b_flag b .|. flagFailsQC }
    p = ent >= r * total
    
    counts = [ count x $ b_seq b | x <- properBases ]
    total = fromIntegral $ V.length $ b_seq b
    ent   = sum [ fromIntegral c * log (total / fromIntegral c) | c <- counts, c /= 0 ] / log 2

-- | Filter on average quality.  Reads without quality string pass.
{-# INLINE qualityAverage #-}
qualityAverage :: Int -> QualFilter
qualityAverage q b = if p then b else b'
  where
    b' = setQualFlag 'Q' $ b { b_flag = b_flag b .|. flagFailsQC }
    p  = let total = S.foldl' (\a x -> a + fromIntegral x) 0 $ b_qual b
         in total >= q * S.length (b_qual b)

-- | Filter on minimum quality.  In @qualityMinimum n q@, a read passes
-- if it has no more than @n@ bases with quality less than @q@.  Reads
-- without quality string pass.
{-# INLINE qualityMinimum #-}
qualityMinimum :: Int -> Word8 -> QualFilter
qualityMinimum n q b = if p then b else b'
  where
    b' = setQualFlag 'Q' $ b { b_flag = b_flag b .|. flagFailsQC }
    p  = S.length (S.filter (< q) (b_qual b)) <= n


-- | Convert quality scores from old Illumina scale (different formula
-- and offset 64 in FastQ).
qualityFromOldIllumina :: BamRec -> BamRec
qualityFromOldIllumina b = b { b_qual = S.map conv $ b_qual b }
  where
    conv :: Word8 -> Word8
    conv s = let s' :: Double
                 s' = exp $ log 10 * (fromIntegral s - 31) / (-10)
                 p  = s' / (1+s')
                 q  = - 10 * log p / log 10
             in round q    
                 
-- | Convert quality scores from new Illumina scale (standard formula
-- but offset 64 in FastQ).
qualityFromNewIllumina :: BamRec -> BamRec
qualityFromNewIllumina b = b { b_qual = S.map (subtract 31) $ b_qual b }


