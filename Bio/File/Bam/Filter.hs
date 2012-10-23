-- | Quality filters adapted from old pipeline.
--
-- TODO: - "SAGE" filter (enforce 17nt reads)?
--       - "rim job" (try to detect glass border)?

{-# LANGUAGE BangPatterns #-}
module Bio.File.Bam.Filter (
    filterPairs, LoneMates(..), QualFilter,
    complexSimple, complexEntropy,
    qualityAverage, qualityMinimum,
    qualityFromOldIllumina, qualityFromNewIllumina
                           ) where

import Bio.File.Bam
import Bio.Iteratee
import Data.Bits
import Data.Word ( Word8 )

import qualified Data.Iteratee    as I
import qualified Data.ByteString  as S

data LoneMates = LoneFail | LoneDrop

-- | A filter applied to pairs of reads.  We supply a predicate to be
-- applied to single reads and one to be applied to pairs.  This
-- function causes an error if a lone mate is hit.  If this is run on a
-- file sorted by coordinate, an error is almost guaranteed.

filterPairs :: Monad m => (BamRec -> Bool) 
                       -> (BamRec -> BamRec -> Bool)
                       -> LoneMates
                       -> Enumeratee [BamRec] [BamRec] m a
filterPairs ps pp lm = eneeCheckIfDone step
  where
    step k = I.tryHead >>= step' k
    step' k Nothing = return $ liftI k
    step' k (Just b)
        | isPaired b = I.tryHead >>= step'' k b
        | otherwise  = if ps b then eneeCheckIfDone step . k $ Chunk [b] else step k

    step'' k b Nothing = case lm of LoneFail -> fail $ "lone mate " ++ show (b_qname b)
                                    LoneDrop -> step k

    step'' k b (Just c)
        | b_rname b /= b_rname c || not (isPaired c) = case lm of LoneFail -> fail $ "lone mate " ++ show (b_qname b)
                                                                  LoneDrop -> step' k (Just c)
        | isFirstMate b && isSecondMate c = step''' k b c
        | isFirstMate c && isSecondMate b = step''' k c b
        | otherwise = case lm of LoneFail -> fail $ "strange pair " ++ show (b_qname b)
                                 LoneDrop -> step' k (Just c)

    step''' k b c = if pp b c then eneeCheckIfDone step . k $ Chunk [b,c] else step k        


-- | A quality filter is simply a transformation on @BamRec@s.  By
-- convention, quality filters should set @flagFailsQC@, a further step
-- can then remove the failed reads.  Filtering of individual reads
-- tends to result in mate pairs with inconsistent flags, which in turn
-- will result in lone mates and all sort of troubles with programs that
-- expect non-broken BAM files.  It is therefore recommended to use
-- @pairFilter@ with suitable predicates to do the post processing.

type QualFilter = BamRec -> BamRec

-- | Simple complexity filter aka "Nancy Filter".  A read is considered
-- not-sufficiently-complex if the most common base accounts for greater
-- than the @cutoff@ fraction of all non-N bases.
{-# INLINE complexSimple #-}
complexSimple :: Double -> QualFilter
complexSimple r b = if p then b else b'
  where
    b' = setQualFlag 'C' $ b { b_flag = b_flag b .|. flagFailsQC }
    p  = let counts = [ length $ filter ((==) x) (b_seq b) | x <- properBases ]
             lim = floor $ r * fromIntegral (sum counts)
         in all (<= lim) counts

-- | Filter on order zero empirical entropy.  Entropy per base must be
-- greater than cutoff.
{-# INLINE complexEntropy #-}
complexEntropy :: Double -> QualFilter
complexEntropy r b = if p then b else b'
  where
    b' = setQualFlag 'C' $ b { b_flag = b_flag b .|. flagFailsQC }
    p = let counts = [ fromIntegral $ length $ filter ((==) x) (b_seq b) | x <- properBases ]
            total = fromIntegral $ length $ b_seq b
            ent   = sum [ c * log (total / c) | c <- counts ] / log 2
        in ent >= r * total

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


