{-# OPTIONS_GHC -funbox-strict-fields #-}
module Bio.File.Pileup where

import Bio.File.Bam.Raw
import Data.Strict.Tuple

-- Pileup of BAM records.  The ingredients for simple genotyping.

type Double4 = Double :!: Double :!: Double :!: Double

-- | The primitive pieces for genotype calling:  A position, a base
-- represented as four probabilities, an inserted sequence, and the
-- length of a deleted sequence.  The logic is that we look at a base
-- followed by some indel, and all those indels are combined into a
-- single insertion and a single deletion.
data PrimChunks 
    = PC { position :: !Int                       -- position on chromosome; we may need to skip forwards
         , probs :: !Double4                      -- four probabilities instead of a call
         , insertion :: [Nucleotide :!: Word8]    -- possible insertion
         , deletion :: !Int                       -- possible deletion
         , more_chunks :: PrimChunks              -- and the remainder
    | EndOfRead


-- | Decomposes a BAM record.  We pick apart the CIGAR field, and
-- combine it with sequence and quality as appropriate.  We ignore the
-- MD field, even if it is present.  Clipped bases are removed,

decompose :: BamRaw -> PrimChunks
decompose br = nextBase (br_pos br) 0 0
  where
    !max_cig = br_cig_len br
    !max_seq = br_l_seq br

    nextBase !pos !is !ic !io | is >= max_seq || ic >= max_cig = EndOfRead
    nextBase !pos !is !ic !io = case co of
        Ins ->             nextBase  pos     (is+cl) (ic+1) 0    -- technically an error, we skip it
        Del ->             nextBase (pos+cl)  is     (ic+1) 0    -- technically an error, we skip it
        Nop ->             nextBase (pos+cl)  is     (ic+1) 0
        SMa ->             nextBase  pos     (is+cl) (ic+1) 0 
        HMa ->             nextBase  pos      is     (ic+1) 0
        Pad ->             nextBase  pos      is     (ic+1) 0
        Mat | io == cl  -> nextBase  pos      is     (ic+1) 0
            | otherwise -> nb $ nextIndel (pos+1) (is+1) ic (io+1) 
      where
        !co = cigar_op  $ br_cigar_at ic br
        !cl = cigar_len $ br_cigar_at ic br

        nb k = let !nc = br_seq_at is br
                   !qu = fromIntegral (br_seq_at is br) + 4.77
             in case () of
                nc == nucA -> k $  0 :!: qu :!: qu :!: qu
                nc == nucC -> k $ qu :!:  0 :!: qu :!: qu
                nc == nucG -> k $ qu :!: qu :!:  0 :!: qu
                nc == nucT -> k $ qu :!: qu :!: qu :!:  0
                otherwise  -> nextBase (pos+1) (is+1) ic (io+1) 

    nextIndel !pos !is !ic !io | is >= max_seq || ic >= max_cig = Invariant EndOfRead
    nextIndel !pos !is !ic !io = case co of -- XXX
        Ins ->             nextBase (is+cl) (ic+1) 0    -- technically an error, we skip it
        Del ->             nextBase  is     (ic+1) 0    -- technically an error, we skip it
        Nop -> Skip cl $   nextIndel is     (ic+1) 0
        SMa ->             nextBase (is+cl) (ic+1) 0 
        HMa ->             nextBase  is     (ic+1) 0
        Pad ->             nextBase  is     (ic+1) 0
        Mat | io == cl  -> nextBase  is     (ic+1) 0
            | otherwise -> nb $ nextIndel (is+1) ic (io+1)
      where
        !co = cigar_op  $ br_cigar_at ic br
        !cl = cigar_len $ br_cigar_at ic br


    let co = cigar_op
        | ld != 0 = MissingBase $ nextIndel is ic (ld-1)
        | br_cigar_at 

    ...





-- nucleotide base along with some sort of error estimate.  
-- | Break a record into pieces useful for genotype calling.
