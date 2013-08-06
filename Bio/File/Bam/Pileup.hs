{-# OPTIONS_GHC -funbox-strict-fields #-}
module Bio.File.Pileup where

import Bio.File.Bam.Raw

-- Pileup of BAM records.  The ingredients for simple genotyping.

-- | The primitive pieces for genotype calling.  We alternate between a
-- base (which may be missing) and and indel (which is often missing).
-- Two co-recursive data types are used for representation, @PrimChunks@
-- and @PrimChunks'@.
data PrimChunks = 
    = SingleBase !Double !Double !Double !Double PrimChunks'
        -- ^ One base, represented as logs probabilities for A, C, G, T
    | Skip !Int PrimChunks'
        -- ^ No data here, because of a deletion or a gap in the CIGAR.
        -- Comes with number of bases to skip.
    | EndOfRead

data PrimChunks'
    = Variant !Int SeqQual PrimChunks
        -- ^ Length of deletion, then inserted sequence with qualities
    | Invariant PrimChunks

data SeqQual = SQ !Nucleotide !Word8 SeqQual | EndOfSeq


-- | Decomposes a BAM record.  We pick apart the CIGAR field, and
-- combine it with sequence and quality as appropriate.  We ignore the
-- MD field, even if it is present.  Clipped bases are removed,

decompose :: BamRaw -> PrimChunks
decompose br = nextBase 0 0
  where
    !max_cig = br_cig_len br
    !max_seq = br_l_seq br

    nextBase !is !ic !io | is >= max_seq || ic >= max_cig = EndOfRead
    nextBase !is !ic !io = case co of
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

        nb = let !nc = br_seq_at is br
                 !qu = fromIntegral (br_seq_at is br) + 4.77
             in case () of
                nc == nucA -> SingleBase 0 qu qu qu 
                nc == nucC -> SingleBase qu 0 qu qu
                nc == nucG -> SingleBase qu qu 0 qu
                nc == nucT -> SingleBase qu qu qu 0
                otherwise  -> Skip 1

    nextIndel !is !ic !io | is >= max_seq || ic >= max_cig = Invariant EndOfRead
    nextIndel !is !ic !io = case co of
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
