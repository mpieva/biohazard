{-# LANGUAGE BangPatterns #-}
{-# OPTIONS_GHC -funbox-strict-fields #-}
module Bio.Bam.Pileup where

-- ^ Genotype Calling:  like Samtools(?), but for aDNA
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
import Bio.Bam.Header
import Bio.Bam.Raw

import Data.List ( intersperse )
import Numeric ( showFFloat )

import qualified Data.ByteString as B

-- | For probabilities for bases @A, C, G, T@.
data Double4 = D4 !Double !Double !Double !Double

instance Show Double4 where
    showsPrec _ (D4 a c g t) = foldr (.) id . intersperse ((:) ' ') $ map (showFFloat (Just 2)) [a,c,g,t]

-- | The primitive pieces for genotype calling:  A position, a base
-- represented as four probabilities, an inserted sequence, and the
-- length of a deleted sequence.  The logic is that we look at a base
-- followed by some indel, and all those indels are combined into a
-- single insertion and a single deletion.
data PrimChunks = Seek !Int PrimBase                            -- ^ skip to position (at start or after N operation)
                | Indel !Int [(Nucleotide, Qual)] PrimBase      -- ^ observed deletion and insertion between two bases
                | EndOfRead                                     -- ^ nothing anymore
  deriving Show

data PrimBase = Base !Double4 PrimChunks                     -- ^ four probabilities instead of a certain base
  deriving Show


-- | Decomposes a BAM record into chunks suitable for piling up.  We
-- pick apart the CIGAR field, and combine it with sequence and quality
-- as appropriate.  We ignore the @MD@ field, even if it is present.
-- Clipped bases are removed/skipped as appropriate.

decompose :: DamageModel -> BamRaw -> PrimChunks
decompose dm br
    | br_isUnmapped br || br_rname br == invalidRefseq = EndOfRead
    | otherwise = firstBase (br_pos br) 0 0
  where
    !max_cig = br_n_cigar_op br
    !max_seq = br_l_seq br
    !mapq    = br_mapq br
    !baq     = br_extAsString "BQ" br

    -- This will compute the effective quality.  As far as I can see
    -- from the BAM spec V1.4, the qualities that matter are QUAL, MAPQ,
    -- and BAQ.  If QUAL is invalid, we replace it (arbitrarily) with
    -- 23 (assuming a rather conservative error rate of ~0.5%), BAQ is
    -- added to QUAL, and MAPQ is an upper limit for effective quality.
    get_seq :: (Nucleotide -> Qual -> a) -> Int -> a
    get_seq f i = f n q''
      where
        !n = br_seq_at br i                                         -- nucleotide
        !q = case br_qual_at br i of 0xff -> 30 ; x -> x            -- quality; invalid (0xff) becomes 30
        !q' | i > B.length baq = q                                  -- no BAQ available
            | otherwise = q + Q (B.index baq i - 64)                -- else correct for BAQ
        !q'' = min q' mapq                                          -- use MAPQ as upper limit

    -- Look for first base following the read's start or a gap (CIGAR
    -- code N).  Indels are skipped, since these are either bugs in the
    -- aligner or the aligner getting rid of essentially unalignable
    -- bases.
    firstBase :: Int -> Int -> Int -> PrimChunks
    firstBase !pos !is !ic
        | is >= max_seq || ic >= max_cig = EndOfRead
        | otherwise = case br_cigar_at br ic of
            (Ins,cl) ->            firstBase  pos (cl+is) (ic+1)
            (Del,cl) ->            firstBase (pos+cl) is  (ic+1)
            (Nop,cl) ->            firstBase (pos+cl) is  (ic+1)
            (SMa,cl) ->            firstBase  pos (cl+is) (ic+1)
            (HMa, _) ->            firstBase  pos     is  (ic+1)
            (Pad, _) ->            firstBase  pos     is  (ic+1)
            (Mat, 0) ->            firstBase  pos     is  (ic+1)
            (Mat, _) -> Seek pos $ nextBase   pos     is   ic 0


    -- Generate probabilities for the next base.  When this gets called,
    -- we are looking at an M CIGAR operation and all the subindices are
    -- valid.
    nextBase :: Int -> Int -> Int -> Int -> PrimBase
    nextBase !pos !is !ic !io = Base (get_seq (dm is) is) $ nextIndel  [] 0 (pos+1) (is+1) ic (io+1)


    -- Look for the next indel after a base.  We collect all indels (I
    -- and D codes) into one combined operation.  If we hit N, we drop
    -- it (indels next to a gap indicate trouble).  Other stuff is
    -- skipped: we could check for stuff that isn't valid in the middle
    -- of a read (H and S), but then what would we do anyway?  Just
    -- ignoring it is much easier and arguably as correct.
    nextIndel :: [[(Nucleotide,Qual)]] -> Int -> Int -> Int -> Int -> Int -> PrimChunks
    nextIndel ins del !pos !is !ic !io
        | is >= max_seq || ic >= max_cig = EndOfRead
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
        isq cl = [ get_seq (,) i | i <- [is..is+cl-1] ] : ins


-- | A 'DamageModel' is a function that gives probabilities for all
-- possible four bases, given the sequenced base, the quality, and the
-- position in the read.  That means it can't take the actual alignment
-- into account... but nobody seems too keen on doing that anyway.
type DamageModel = Int              -- ^ position in read
                -> Nucleotide       -- ^ base
                -> Qual             -- ^ quality score
                -> Double4          -- ^ results in four probabilities

-- | 'DamageModel' for undamaged DNA.  The probabilities follow directly
-- from the quality score.  This needs elaboration to see what to do
-- with amibiguity codes.
noDamage :: DamageModel
noDamage _ b (Q q) | b == nucA = D4 0 p p p
                   | b == nucC = D4 p 0 p p
                   | b == nucG = D4 p p 0 p
                   | b == nucT = D4 p p p 0
                   | otherwise = D4 0 0 0 0
  where !p = fromIntegral q + 4.77


-- | A variant call consists of a position, some measure of qualities,
-- and the actual call of either a SNV or an Indel.  SNV and Indel calls
-- will alternate if both are available, but enforcing this is too
-- cumbersome.  Also, most Indel calls won't be generated for lack of
-- evidence anyway.

data VarCall = VarCall { vc_refseq     :: !Refseq
                       , vc_pos        :: !Int
                       , vc_depth      :: !Int         -- number of contributing reads
                       , vc_mapq0      :: !Int         -- number of contributing reads with MAPQ==0
                       , vc_sum_mapq   :: !Int         -- sum of map qualities of contributring reads
                       , vc_sum_mapq2  :: !Int         -- sum of squared map qualities of contributing reads
                       , vc_call       :: Varcall' }   -- variant call, whatever that means

data VarCall' = SnvCall   { pl   :: !UVector Int Qual }     -- PL values in dB
              | IndelCall { pl   :: !UVector Int Qual       -- PL values in dB
                          , vars :: [(Int,[Nucleotide])]    -- variant: number of deletions, inserted sequence)


-- | Running pileup results in a series of piles.  A 'Pile' has a
-- position, and contents depending on whether bases or indels were
-- piled up.

data Pile = Pile { p_refseq :: !Refseq
                 , p_pos    :: !Int
                 , p_pile   :: Either BasePile IndelPile }

-- | A 'BasePile' contains map qualities for the contributing reads and
-- four probabilities for the four possible bases.

type BasePile = [ (Qual, Double4) ]

-- | An 'IndelPile' contains mapq qualities for the contributing reads,
-- the amount of sequence deleted for each, and the sequence inserted
-- for each.  Note that a variant can both delete and insert (if the
-- aligner thinks reporting the changes in this way makes sense).

type IndelPile = [ (Qual, Int, [Nucleotide]) ]


-- | A gentoype caller receives a list of bases at a column and is
-- supposed to turn them into something... PL values come to mind.  XXX
-- Trouble is, the 'something' may well change with the concrete caller.
-- (This used to have a random component.  I'm not going to support
-- sillyness like that anymore:  *this* programm will not throw out
-- information.  Use another one if that's your cup of tea.)

-- type Conscall = [Double4] -> Double10

-- | The simplest genotype caller imaginable:  majority call, weighted by quality scores.
{-majorityCall :: Conscall
majorityCall ns gen = case qns of
    (n1,q1) : (_,q2) : _ -> Just n1 :!: q1-q2 :!: gen
    (n1,q1) : []         -> Just n1 :!:    q1 :!: gen
    []                   -> Nothing :!:     0 :!: gen
  where
    qns = sortBy (\(_,a) (_,b) -> b `compare` a) $ M.toList $
          foldl' (\m (n:!:q) -> M.insertWith' (+) n q m) M.empty ns -}


-- | The pileup enumeratee takes 'BamRec's, decomposes them, interleaves
-- the pieces appropriately, and generates 'Pile's.  The output will
-- contain at most one 'BasePile' and one 'IndelPile' for each position,
-- piles are sorted by position.

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
            !gs'  = maybe gs (: gs) mbq -}

-- Applies a consensus caller.  Threads an RNG through, counts length.
-- Useless, except for reference.
{- appConscall :: Monad m => StdGen -> Conscall -> Enumeratee [Pile] [Column] m a
appConscall gen0 ccall = eneeCheckIfDone (liftI . go gen0)
  where
    go  _  out (EOF      mx) = idone (liftI out) $ EOF mx
    go gen out (Chunk piles) = eneeCheckIfDone (liftI . go gen') . out . Chunk $ cols
      where (gen', cols) = mapAccumL step gen piles

    step gen (Pile s p b grs) = case ccall grs gen of
        mn :!: q :!: gen' -> let !l = length grs in (gen', Column s p [Just b :!: 30 :!: 1, mn :!: q :!: l])
-}
