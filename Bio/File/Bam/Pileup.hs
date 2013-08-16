{-# LANGUAGE BangPatterns #-}
{-# OPTIONS_GHC -funbox-strict-fields #-}
module Bio.File.Bam.Pileup where

import Bio.Base
import Bio.File.Bam
import Bio.File.Bam.Raw

import Data.List ( intersperse )
import Numeric ( showFFloat )

-- Pileup of BAM records.  The ingredients for simple genotyping.

data Double4 = D4 !Double !Double !Double !Double

instance Show Double4 where
    showsPrec _ (D4 a c g t) = foldr (.) id . intersperse ((:) ' ') $ map (showFFloat (Just 2)) [a,c,g,t]

-- | The primitive pieces for genotype calling:  A position, a base
-- represented as four probabilities, an inserted sequence, and the
-- length of a deleted sequence.  The logic is that we look at a base
-- followed by some indel, and all those indels are combined into a
-- single insertion and a single deletion.
data PrimChunks = Seek !Int PrimBase                            -- skip to position (at start or after N operation)
                | Indel !Int [(Nucleotide, Int)] PrimBase       -- observed deletion and insertion
                | EndOfRead                                     -- nothing anymore
  deriving Show

data PrimBase = Base !Double4 PrimChunks                     -- four probabilities instead of a base
  deriving Show


-- | Decomposes a BAM record.  We pick apart the CIGAR field, and
-- combine it with sequence and quality as appropriate.  We ignore the
-- MD field, even if it is present.  Clipped bases are removed,

decompose :: DamageModel -> BamRaw -> PrimChunks
decompose dm br = firstBase (br_pos br) 0 0
  where
    !max_cig = br_n_cigar_op br
    !max_seq = br_l_seq br

    firstBase :: Int -> Int -> Int -> PrimChunks
    firstBase !pos !is !ic
        | is >= max_seq || ic >= max_cig = EndOfRead
        | otherwise = case br_cigar_at br ic of
            (Ins,cl) ->            firstBase  pos     (is+cl) (ic+1)    -- indel at beginning is technically an error
            (Del,cl) ->            firstBase (pos+cl)  is     (ic+1)    -- indel at beginning is technically an error
            (Nop,cl) ->            firstBase (pos+cl)  is     (ic+1)
            (SMa,cl) ->            firstBase  pos     (is+cl) (ic+1)
            (HMa, _) ->            firstBase  pos      is     (ic+1)
            (Pad, _) ->            firstBase  pos      is     (ic+1)
            (Mat, 0) ->            firstBase  pos      is     (ic+1)
            (Mat, _) -> Seek pos $ nextBase   pos      is      ic    0


    nextBase :: Int -> Int -> Int -> Int -> PrimBase
    nextBase !pos !is !ic !io = Base bq $ nextIndel  [] 0 (pos+1) (is+1) ic (io+1)
      where
        !bq = dm is (br_seq_at br is) (br_qual_at br is)


    nextIndel :: [[(Nucleotide,Int)]] -> Int -> Int -> Int -> Int -> Int -> PrimChunks
    nextIndel ins del !pos !is !ic !io
        | is >= max_seq || ic >= max_cig = EndOfRead                    -- indel at end is technically an error
        | otherwise = case br_cigar_at br ic of
            (Ins,cl) ->             nextIndel (isq cl) del   pos (cl+is) (ic+1) 0
            (Del,cl) ->             nextIndel  ins (cl+del) (pos+cl) is  (ic+1) 0
            (SMa,cl) ->             nextIndel  ins     del   pos (cl+is) (ic+1) 0
            (Pad, _) ->             nextIndel  ins     del   pos     is  (ic+1) 0
            (HMa, _) ->             nextIndel  ins     del   pos     is  (ic+1) 0
            (Mat,cl) | io == cl  -> nextIndel  ins     del   pos     is  (ic+1) 0
                     | otherwise -> Indel del out $ nextBase pos     is   ic   io   -- ends up generating a 'Base'
            (Nop,cl) ->             firstBase               (pos+cl) is  (ic+1)     -- ends up generating a 'Seek'
      where
        out    = concat $ reverse ins
        isq cl = [ (br_seq_at br i, br_qual_at br i) | i <- [is..is+cl-1] ] : ins


type DamageModel = Int              -- ^ position in read
                -> Nucleotide       -- ^ base
                -> Int              -- ^ quality score
                -> Double4          -- ^ results in four probabilities

-- | 'DamageModel' for undamaged DNA.  The probabilities follow directly
-- from the quality score.  This needs elaboration to see what to do
-- with amibiguity codes.
noDamage :: DamageModel
noDamage _ b q | b == nucA = D4 0 p p p
               | b == nucC = D4 p 0 p p
               | b == nucG = D4 p p 0 p
               | b == nucT = D4 p p p 0
               | otherwise = D4 0 0 0 0
  where !p = fromIntegral q + 4.77


-- pileup :: Monad m => Enumeratee [PrimChunks] [Pile] m a
