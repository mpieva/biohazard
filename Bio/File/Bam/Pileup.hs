{-# LANGUAGE BangPatterns #-}
{-# OPTIONS_GHC -funbox-strict-fields #-}
module Bio.File.Bam.Pileup where

import Bio.Base
import Bio.File.Bam
import Bio.File.Bam.Raw

-- Pileup of BAM records.  The ingredients for simple genotyping.

data Double4 = D4 !Double !Double !Double !Double

-- | The primitive pieces for genotype calling:  A position, a base
-- represented as four probabilities, an inserted sequence, and the
-- length of a deleted sequence.  The logic is that we look at a base
-- followed by some indel, and all those indels are combined into a
-- single insertion and a single deletion.
data PrimChunks 
    = PC { position    :: !Int                      -- position on chromosome; we may need to skip forwards
         , probs       :: !Double4                  -- four probabilities instead of a call
         , insertion   :: [(Nucleotide, Int)]       -- possible insertion
         , deletion    :: !Int                      -- possible deletion
         , more_chunks :: PrimChunks }              -- and the remainder
    | EndOfRead


-- | Decomposes a BAM record.  We pick apart the CIGAR field, and
-- combine it with sequence and quality as appropriate.  We ignore the
-- MD field, even if it is present.  Clipped bases are removed,

decompose :: BamRaw -> PrimChunks
decompose br = nextBase (br_pos br) 0 0 0
  where
    !max_cig = br_n_cigar_op br
    !max_seq = br_l_seq br

    nextBase :: Int -> Int -> Int -> Int -> PrimChunks
    nextBase !pos !is !ic !io
        | is >= max_seq || ic >= max_cig = EndOfRead
        | otherwise = case co of
            Ins ->             nextBase  pos     (is+cl) (ic+1) 0    -- technically an error, we skip it
            Del ->             nextBase (pos+cl)  is     (ic+1) 0    -- technically an error, we skip it
            Nop ->             nextBase (pos+cl)  is     (ic+1) 0
            SMa ->             nextBase  pos     (is+cl) (ic+1) 0 
            HMa ->             nextBase  pos      is     (ic+1) 0
            Pad ->             nextBase  pos      is     (ic+1) 0
            Mat | io == cl  -> nextBase  pos      is     (ic+1) 0
                | otherwise -> nextIndel (PC pos bq) [] 0 (pos+1) (is+1) ic (io+1) 
      where
        !co = cigar_op  $ br_cigar_at br ic
        !cl = cigar_len $ br_cigar_at br ic

        bq | nc == nucA = D4  0 qu qu qu
           | nc == nucC = D4 qu  0 qu qu
           | nc == nucG = D4 qu qu  0 qu
           | nc == nucT = D4 qu qu qu  0
           | otherwise  = D4  0  0  0  0       -- this needs elaboration...
          where
            !nc = br_seq_at br is
            !qu = fromIntegral (br_qual_at br is) + 4.77
               

    nextIndel :: ([(Nucleotide, Int)] -> Int -> PrimChunks -> PrimChunks) 
              -> [[(Nucleotide, Int)]] -> Int -> Int -> Int -> Int -> Int -> PrimChunks
    nextIndel k ins del !pos !is !ic !io
        | is >= max_seq || ic >= max_cig = k out del EndOfRead
        | otherwise = case co of
            Ins -> let isq = [ (br_seq_at br i, br_qual_at br i) | i <- [is..is+cl-1] ]
                   in nextIndel k (isq:ins) del      pos     (is+cl) (ic+1) 0
            Del ->    nextIndel k      ins (del+cl) (pos+cl)  is     (ic+1) 0
            Pad ->    nextIndel k      ins  del      pos      is     (ic+1) 0
            Nop -> k out del $ nextBase (pos+cl) is     (ic+1) 0
            SMa -> k out del $ nextBase  pos    (is+cl) (ic+1) 0
            HMa -> k out del $ nextBase  pos     is     (ic+1) 0 
            Mat | io == cl  -> nextIndel k ins del            pos is (ic+1) 0
                | otherwise ->           k out del $ nextBase pos is  ic   io
      where
        !co = cigar_op  $ br_cigar_at br ic
        !cl = cigar_len $ br_cigar_at br ic
        out = concat $ reverse ins


-- nucleotide base along with some sort of error estimate.  
-- | Break a record into pieces useful for genotype calling.
