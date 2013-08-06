{-# LANGUAGE BangPatterns, PatternGuards #-}
module Bio.Pileup where

-- Genotype Calling:  for aDNA, but otherwise rather simple minded
--
-- The goal for this module is to call haploid and diploid single
-- nucleotide variants the best way we can, including support for aDNA.
-- Indel calling is outr of scope, we only do it "on the side".
--
-- The cleanest way to call genotypes under all circumstances is
-- probably the Dindel approach:  define candidate haplotypes, align
-- each read to each haplotype, then call the likely haplotypes with a
-- quality derived from the quality scores.  This approach neatly
-- integrates indel calling with ancient DNA and makes a separate indel
-- realigner redundant.  However, it's rather expensive in that it
-- requires inclusion of an aligner, and we'd need an aligner that is
-- compatible with the chosen error model, which can be hard.
--
-- Here we'll take a short cut:  We do not really call indels; these
-- variants are collected and are assigned a affine score.  This works
-- best if indels are 'left-aligned' first.  In theory, one indel
-- variant could be another indel variant with a sequencing error---we
-- ignore that.  Once indels are taken care off, SNVs are treated
-- separately as independent columns of the pileup.
--
-- For aDNA, we need a substitution probability.  We have three options:
-- use an empirically determined PSSM, use an arithmetically defined
-- PSSM based on the Johnson model, use a context sensitive PSSM based
-- on the Johnson model and an alignment.  Using Dindel, actual
-- substitutions relative to a called haplotype would be taken into
-- account.  Since we're not going to do that, taking alignments into
-- account is difficult, somewhat approximatem, and not worth the
-- hassle.
--
-- Regarding the error model, there's a choice between samtools/maq or
-- the naive model everybody else uses.  Naive is easy to marry to aDNA,
-- samtools is (probably) better.  Either way, we introduce a number of
-- parameters (eta and kappa for samtools, lambda and p for Johnson).
-- Running a maximum likehood fit for those may be valuable.  It would
-- be cool, if we could do that without rerunning the complete genotype
-- caller, but it's not a priority.
--
-- So, outline of the genotype caller:  We read BAM (not filtering at
-- all, that's somebody else's problem, but we might want to split by
-- read group).  We will scan each read's CIGAR line in concert with the
-- sequence and effective quality.  Effective quality is the lowest
-- available quality score of QUAL, MAPQ, and BQ.  For aDNA calling, the
-- base is transformed into four likelihoods based on the aDNA
-- substitution matrix.
--
-- So, either way, we need something like "pileup", where indel variants
-- are collected as they are (any length), while matches are piled up.


import Bio.Base
import Bio.File.Bam
import Bio.Iteratee
import Control.Monad
import Data.List
import Data.Sequence ((|>))
import System.Random

import qualified Data.ByteString        as S
import qualified Data.Iteratee.ListLike as L
import qualified Data.Foldable          as Z
import qualified Data.Map               as M
import qualified Data.Sequence          as Z
import qualified Data.Vector.Unboxed    as V


data VarCall 
    = SnvCall { pl :: [ Int ] }             -- PL values in Phred scale, 4 haploid, 10 diploid
    | IndelCall { vars :: [[Nucleotide]]
                , pl :: [ Int ] }           -- PL values, number depends on number of alleles
    

-- Read BAM, pile up the alignments so we get a collection of bases per
-- column.  We cook every read into a list of fragments:

data Frags = Wait Position Frags        -- not at the right position
           | Insert Sequence Frags'     -- inserted sequence
           | Delete Sequence Frags'     -- deleted sequence
           | NoIndel Frags'             -- no indel here
           | EndOfRead

data Frags' = SkipBase Frags                    -- was deleted
            | MatchBase Nucleotide Int Frags    -- base and quality

-- | Cook a read into a list of fragments.  Fragments alternate between
-- indels and bases, bases are annotated with an effective quality
-- score.  We skip over reads that aren't aligned.
cookRead :: BamRaw -> Frags
cookRead br | br_isUnmapped br || br_rname br == invalidRefseq = EndOfRead
            | otherwise = Wait (Pos (br_rname br) (br_pos br)) frags
  where
    baq = br_extAsString "BQ" br
    eff_qual | S.null baq = map (min $ br_mapq br) (br_qual br)
             | otherwise  = map (min $ br_mapq br) $ zipWith
                                (\q o -> q - (o - 64))
                                (br_qual br) baq

    frags = go (br_seq br) (br_qual br) (



-- Consensus callers.
-- Receive a list of bases with qualities and a random generator.
-- Return a quality and a new RNG, might return a base.
type Conscall = [Nucleotide :!: Int] -> StdGen -> Maybe Nucleotide :!: Int :!: StdGen

-- Majority call, weighted by quality scores.
majorityCall :: Conscall
majorityCall ns gen = case qns of
    (n1,q1) : (_,q2) : _ -> Just n1 :!: q1-q2 :!: gen
    (n1,q1) : []         -> Just n1 :!:    q1 :!: gen
    []                   -> Nothing :!:     0 :!: gen
  where
    qns = sortBy (\(_,a) (_,b) -> b `compare` a) $ M.toList $ 
          foldl' (\m (n:!:q) -> M.insertWith' (+) n q m) M.empty ns

-- XXX missing: diploid consensus caller.  Must call SNPs to be able to
-- ask for allele skew.  Means I need an encoding for ambiguity codes.



-- Run consensus caller, but also filter for sensible coverage
withinCoverage :: (Int,Int) -> Conscall -> Conscall
withinCoverage (low,high) call ns gen 
    | length ns < low || length ns > high = Nothing :!: 0 :!: gen
    | otherwise                           = call ns gen


-- Pileup overlapping reads _on_a_single_chromosome_.  Produces a stream
-- of positions with reference base and list of observations.  Ignores
-- gaps; observed gaps don't appear in the output, inserted bases are
-- dropped.
-- (Note: we can't easily get reference data outside an alignment.  *If*
-- we want to filter for CpG sites, we may need to drop the first/last
-- base unconditionally.)

data BRead = BRead { r_name     :: !Seqid
                   , r_sequence :: [Nucleotide]
                   , r_quality  :: [Int]
                   , r_cigar    :: [(CigOp,Int)]
                   , r_md       :: Maybe [MdOp] }
    deriving Show

-- | Unpacks reads.  Reads must have MD fields, else they are silently
-- dropped.  Missing quality scores are assumed to be 30, quality scores
-- are limited to MAPQ.
unpackRead :: BamRec -> BRead
unpackRead b = BRead (b_qname b) (V.toList $ b_seq b) qs (unCigar $ b_cigar b) (P.maybe Nothing Just $ getMd b)
   where qs | S.null (b_qual b) = replicate (V.length $ b_seq b) (min 30 (b_mapq b))
            | otherwise         = map (min (b_mapq b) . fromIntegral) $ S.unpack $ b_qual b

pileup :: Monad m => Enumeratee [BamRec] [Pile] m a
pileup = eneeCheckIfDone (liftI . pileup' invalidRefseq 0 Z.empty)

pileup' :: Monad m 
        => Refseq -> Int -> Z.Seq BRead
        -> (Stream [Pile] -> Iteratee [Pile] m a)
        -> Stream [BamRec] -> Iteratee [BamRec] m (Iteratee [Pile] m a)
pileup' !_ !_ !pile !out (EOF       mx) | Z.null pile  = idone (liftI out) $ EOF mx
pileup' !s !p !pile !out (EOF       mx)                = emit s p pile out $ EOF mx
pileup' !s !p !pile !out (Chunk     [])                = liftI $ pileup' s p pile out
pileup' !s !p !pile !out (Chunk (r:rs))
    | Z.null pile && good                           = pileup' (b_rname r) (b_pos r) (Z.singleton $ unpackRead r) out (Chunk rs)
    | Z.null pile                                   = pileup' s p Z.empty out (Chunk rs)
    | s /= b_rname r || p < b_pos r                 = emit    s p pile out (Chunk (r:rs))
    | good                                          = pileup' s p (pile |> unpackRead r) out (Chunk rs)
    | otherwise                                     = pileup' s p pile out (Chunk rs)
  where
    good = not (isUnmapped r || isFailsQC r)

emit :: Monad m 
     => Refseq -> Int -> Z.Seq BRead
     -> (Stream [Pile] -> Iteratee [Pile] m a)
     -> Stream [BamRec] -> Iteratee [BamRec] m (Iteratee [Pile] m a)
emit !s !p !rs out stream = case m'rb of 
        Just rb -> eneeCheckIfDone cont (out $ Chunk [Pile s p rb grs])
        Nothing -> pileup' s (p+1) rs' out stream
  where
    cont it = pileup' s (p+1) rs' it stream
    m'rb :!: grs :!: rs' = Z.foldl' (\acc rd -> step acc $ unconsRead rd) (Nothing :!: [] :!: Z.empty) rs

    step (mrb :!: gs :!: rs0)                  Nothing  = mrb  :!: gs  :!: rs0
    step (mrb :!: gs :!: rs0) (Just (mb :!: mbq :!: r)) = mrb' :!: gs' :!: rs0 |> r
      where !mrb' = maybe mb Just mrb
            !gs'  = maybe gs (: gs) mbq
    

-- Split off first base at next available reference position(!), return
-- the reference base at that position, the first base & quality, and
-- the shortened read.  Deals with CIGAR, SEQ, QUAL, POS and MD; we
-- don't need the other fields.  If the read turns out to be empty (all
-- CIGAR and MD ops are of length 0, no sequence and quality), @Nothing@
-- is returned.
unconsRead :: BRead -> Maybe (Maybe Nucleotide :!: Maybe (Nucleotide :!: Int) :!: BRead)
unconsRead rd = go (r_cigar rd) (r_sequence rd) (r_quality rd) (r_md rd)
  where
    go ((  _,0):cig)  sq qs md = go cig sq qs md
    go cig sq qs (Just (MdNum  0:md)) = go cig sq qs (Just md)
    go cig sq qs (Just (MdDel []:md)) = go cig sq qs (Just md)

    go ((Ins,n):cig) sq qs md = go cig (drop n sq) (drop n qs) md
    go ((SMa,n):cig) sq qs md = go cig (drop n sq) (drop n qs) md
    go ((Nop,n):cig) sq qs md = r0 (n-1) cig sq qs md
    go ((HMa,_):cig) sq qs md = go cig sq qs md
    go ((Pad,_):cig) sq qs md = go cig sq qs md

    -- deleted reference base aka gap
    go ((Del,n):cig)    sq     qs  (Just (MdDel (y:dd):md)) = r   y  gap 0 ((Del,n-1):cig) sq qs (Just $ MdDel  dd  :md)
    go ((Del,n):cig)    sq     qs                   Nothing = r nucN gap 0 ((Del,n-1):cig) sq qs Nothing
    -- matched, mismatched reference base
    go ((Mat,n):cig) (x:sq) (q:qs) (Just (MdNum   m   :md)) = r   x   x  q ((Mat,n-1):cig) sq qs (Just $ MdNum (m-1):md)
    go ((Mat,n):cig) (x:sq) (q:qs) (Just (MdRep   y   :md)) = r   y   x  q ((Mat,n-1):cig) sq qs (Just md)
    go ((Mat,n):cig) (x:sq) (q:qs)                  Nothing = r nucN  x  q ((Mat,n-1):cig) sq qs Nothing

    -- we might actually decompose the whole read and end up with all empty lists
    go [] [] [] (Just []) = Nothing
    go [] [] []  Nothing  = Nothing

    -- But if only some lists are empty or the operations don't match
    -- up, then we complain.
    go _ _ _ _ = error $ "CIGAR and SEQ inconsistent with each other in " ++ show rd

    r  !n !x !q cig sq !qs md = Just (Just  n :!: Just (x :!: q) :!: BRead (r_name rd) sq qs cig md)
    r0 !n'      cig sq !qs md = Just (Nothing :!:        Nothing :!: BRead (r_name rd) sq qs ((Nop,n'):cig) md)

appConscall :: Monad m => StdGen -> Conscall -> Enumeratee [Pile] [Column] m a
appConscall gen0 ccall = eneeCheckIfDone (liftI . go gen0)
  where
    go  _  out (EOF      mx) = idone (liftI out) $ EOF mx
    go gen out (Chunk piles) = eneeCheckIfDone (liftI . go gen') . out . Chunk $ cols
      where (gen', cols) = mapAccumL step gen piles

    step gen (Pile s p b grs) = case ccall grs gen of
        mn :!: q :!: gen' -> let !l = length grs in (gen', Column s p [Just b :!: 30 :!: 1, mn :!: q :!: l]) 
        

test0 :: IO ()
test0 = getStdGen >>= \gen ->
        enumFile defaultBufSize "foo.bam" >=> run $
        joinI $ decompressBgzf $
        joinI $ decodeBam $ \_ -> 
        joinI $ L.mapStream decodeBamEntry $
        joinI $ pileup $
        joinI $ appConscall gen majorityCall $
        L.mapM_ print


