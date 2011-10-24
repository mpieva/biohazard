{-# LANGUAGE BangPatterns #-}
module Bio.File.Bam.Filter (
    QualFilter, qualityFilterWith, qualityFilterWith',
    complexSimple, complexEntropy,
    qualityAverage, qualityMinimum
                           ) where

-- Quality filters adapted from old pipeline.
-- TODO: - "SAGE" filter (enforce 17nt reads)?
--       - "rim job" (try to detect glass border)?
--       - quality conversion (old Solexa to Phred scale)?

import Bio.File.Bam
import Bio.Iteratee
import Control.Monad.Trans.Class
import Data.Bits
import Data.Word ( Word8 )

import qualified Data.Iteratee    as I
import qualified Data.ByteString  as S

-- | Generic quality filter.  In @qualityFilterWith (p,f)@, a read @r@
-- passes the filter if @p r@ is @True@.  If all reads in a group of
-- consecutive records with the same QNAME pass the filter, the whole
-- group is passed through.  Else @f@ is applied to each read in the
-- group(!) before passing it on and each read in the group(!) is
-- considered to have failed the filter.
-- The return value is @(u,v,it)@, where @u@ is the number of reads that
-- were not marked @isFailsQC@ before and passed the filter, @v@ is the
-- number of reads that were not marked @isFailsQC@ before, and @it@ is
-- the transformed output @Iteratee@.

qualityFilterWith :: Monad m => QualFilter -> Iteratee [BamRec] m a
                  -> Iteratee [BamRec] m (Int, Int, Iteratee [BamRec] m a)
qualityFilterWith (p,f) = joinI . I.groupBy same_qname . liftI . stepQF p f 0 0
  where same_qname a b = b_qname a == b_qname b

-- | Same as @qualityFilterWith@, but reads are not grouped before
-- testing.  This will result in mate pairs with inconsistent flags,
-- which in turn will result in lone mates and all sort of troubles with
-- programs that expect non-broken BAM files.  Present merely for
-- completeness, not because it's useful.
{-# DEPRECATED qualityFilterWith' "You are strongly urged to consider qualityFilterWith instead" #-}
qualityFilterWith' :: Monad m => QualFilter -> Iteratee [BamRec] m a
                   -> Iteratee [BamRec] m (Int, Int, Iteratee [BamRec] m a)
qualityFilterWith' (p,f) = joinI . I.mapStream (:[]) . liftI . stepQF p f 0 0

stepQF :: (Monad m) => (BamRec -> Bool) -> (BamRec -> BamRec)
     -> Int -> Int -> Iteratee [BamRec] m a
     -> Stream [[BamRec]] -> Iteratee [[BamRec]] m (Int, Int, Iteratee [BamRec] m a)
stepQF p f = step    
  where
    step !u !v it (EOF       mx) = idone (u,v,it) (EOF mx)
    step !u !v it (Chunk [    ]) = liftI $ step u v it
    step !u !v it (Chunk (r:rs)) 
        | all p r   = lift (enumPure1Chunk r it) >>= \it' ->
                      step (u+c) (v+c) it' (Chunk rs)
        | otherwise = lift (enumPure1Chunk (map f r) it) >>= \it' ->
                      step u (v+c) it' (Chunk rs)
      where !c = length $ filter (not . isFailsQC) r

type QualFilter = (BamRec->Bool, BamRec->BamRec)

-- | Simple complexity filter aka "Nancy Filter".  A read is considered
-- not-sufficiently-complex if the most common base accounts for greater
-- than the @cutoff@ fraction of all non-N bases.
complexSimple :: Double -> QualFilter
complexSimple r = (p,f)
  where
    f b = b { b_flag = b_flag b .|. flagFailsQC .|. flagLowComplexity }
    p b = let counts = [ length $ filter ((==) x) (b_seq b) | x <- [A,C,G,T] ]
              lim = floor $ r * fromIntegral (sum counts)
          in all (<= lim) counts

-- | Filter on order zero empirical entropy.  Average entropy must be
-- greater than cutoff.
complexEntropy :: Double -> QualFilter
complexEntropy r = (p,f)
  where
    f b = b { b_flag = b_flag b .|. flagFailsQC .|. flagLowComplexity }
    p b = let counts = [ fromIntegral $ length $ filter ((==) x) (b_seq b) | x <- [A,C,G,T] ]
              total = fromIntegral $ length $ b_seq b
              ent   = sum [ c * log (total / c) | c <- counts ] / log 2
          in ent >= r * total

-- | Filter on average quality.  Reads without quality string pass.
qualityAverage :: Int -> QualFilter
qualityAverage q = (p,f)
  where
    f b = b { b_flag = b_flag b .|. flagFailsQC .|. flagLowQuality }
    p b = let total = S.foldl' (\a x -> a + fromIntegral x) 0 $ b_qual b
          in total >= q * S.length (b_qual b)

-- | Filter on minimum quality.  In @qualityMinimum n q@, a read passes
-- if it has no more than @n@ bases with quality less than @q@.  Reads
-- without quality string pass.
qualityMinimum :: Int -> Word8 -> QualFilter
qualityMinimum n q = (p,f)
  where
    f b = b { b_flag = b_flag b .|. flagFailsQC .|. flagLowQuality }
    p b = S.length (S.filter (< q) (b_qual b)) <= n
