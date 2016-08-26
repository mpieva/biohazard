{-# LANGUAGE BangPatterns, Rank2Types, RecordWildCards, OverloadedStrings #-}
{-# OPTIONS_GHC -funbox-strict-fields #-}
module Bio.Bam.Pileup where

import Bio.Adna ( Mat44D )
import Bio.Bam.Header
import Bio.Bam.Rec
import Bio.Iteratee
import Bio.Prelude

import qualified Data.ByteString        as B
import qualified Data.Vector.Generic    as V
import qualified Data.Vector.Unboxed    as U

-- ^ Genotype Calling:  like Samtools(?), but for aDNA
--
-- The goal for this module is to call haploid and diploid single
-- nucleotide variants the best way we can, including support for aDNA.
-- Indel calling is out of scope, we only do it "on the side".
--
-- The cleanest way to call genotypes under all circumstances is
-- probably the /Dindel/ approach:  define candidate haplotypes, align
-- each read to each haplotype, then call the likely haplotypes with a
-- quality derived from the quality scores.  This approach neatly
-- integrates indel calling with ancient DNA and makes a separate indel
-- realigner redundant.  However, it's rather expensive in that it
-- requires inclusion of an aligner, and we'd need an aligner that is
-- compatible with the chosen error model, which might be hard.
--
-- Here we'll take a short cut:  We do not really call indels.  Instead,
-- these variants are collected and are assigned an affine score.  This
-- works best if indels are 'left-aligned' first.  In theory, one indel
-- variant could be another indel variant with a sequencing error---we
-- ignore that possibility for the most part.  Once indels are taken
-- care off, SNVs are treated separately as independent columns of the
-- pileup.
--
-- Regarding the error model, there's a choice between /samtools/ or the
-- naive model everybody else (GATK, Rasmus Nielsen, etc.) uses.  Naive
-- is easy to marry to aDNA, samtools is (probably) better.  Either way,
-- we introduce a number of parameters (@eta@ and @kappa@ for
-- /samtools/, @lambda@, @delta@, @delta_ss@ for /Johnson/).  Running a
-- maximum likehood fit for those may be valuable.  It would be cool, if
-- we could do that without rerunning the complete genotype caller, but
-- it's not a priority.
--
-- So, outline of the genotype caller:  We read BAM (minimally
-- filtering; general filtering is somebody else's problem, but we might
-- want to split by read group).  We will scan each read's CIGAR line in
-- concert with the sequence and effective quality.  Effective quality
-- is the lowest available quality score of QUAL, MAPQ, and BQ.  For
-- aDNA calling, the base is transformed into four likelihoods based on
-- the aDNA substitution matrix.
--
-- So, either way, we need something like "pileup", where indel variants
-- are collected as they are (any length), while matches are piled up.
--
-- Regarding output, we certainly don't want to write VCF or BCF.  (No
-- VCF because it's ugly, no BCF, because the tool support is
-- non-existent.)  It will definitely be something binary.  For the GL
-- values, small floating point formats may make sense: half-precision
-- floating point's representable range would be 6.1E-5 to 6.5E+5, 0.4.4
-- minifloat from Bio.Util goes from 0 to 63488.


-- *TODO*
--
-- * A whole lot of testing.
-- * ML fitting and evaluation of parameters for different possible
--   error and damage models.
-- * Maybe specialize to ploidy one and two.

-- | The primitive pieces for genotype calling:  A position, a base
-- represented as four likelihoods, an inserted sequence, and the
-- length of a deleted sequence.  The logic is that we look at a base
-- followed by some indel, and all those indels are combined into a
-- single insertion and a single deletion.
data PrimChunks = Seek Int PrimBase                             -- ^ skip to position (at start or after N operation)
                | Indel [Nucleotides] [DamagedBase] PrimBase    -- ^ observed deletion and insertion between two bases
                | EndOfRead                                     -- ^ nothing anymore
  deriving Show

data PrimBase = Base { _pb_wait   :: {-# UNPACK #-} !Int        -- ^ number of bases to wait due to a deletion
                     , _pb_likes  :: DamagedBase                -- ^ four likelihoods
                     , _pb_mapq   :: {-# UNPACK #-} !Qual       -- ^ map quality
                     , _pb_rev    :: Bool                       -- ^ reverse strand?
                     , _pb_chunks :: PrimChunks }               -- ^ more chunks
  deriving Show

type PosPrimChunks = (Refseq, Int, PrimChunks)

-- | Represents our knowledge about a certain base, which consists of
-- the base itself (A,C,G,T, encoded as 0..3; no Ns), the quality score
-- (anything that isn't A,C,G,T becomes A with quality 0), and a
-- substitution matrix representing post-mortem but pre-sequencing
-- substitutions.
--
-- Unfortunately, none of this can be rolled into something more simple,
-- because damage and sequencing error behave so differently.

data DamagedBase = DB { db_call :: {-# UNPACK #-} !Nucleotide           -- ^ called base
                      , db_qual :: {-# UNPACK #-} !Qual                 -- ^ quality of called base
                      , db_ref  :: {-# UNPACK #-} !Nucleotides          -- ^ reference base from MD field
                      , db_dmg  :: {-# UNPACK #-} !Mat44D }             -- ^ damage matrix

instance Show DamagedBase where
    showsPrec _ (DB n q r _)
        | nucToNucs n == r = shows n .                     (:) '@' . shows q
        | otherwise        = shows n . (:) '/' . shows r . (:) '@' . shows q


-- | Decomposes a BAM record into chunks suitable for piling up.  We
-- pick apart the CIGAR and MD fields, and combine them with sequence
-- and quality as appropriate.  Clipped bases are removed/skipped as
-- appropriate.  We also do apply a substitution matrix to each base,
-- which must be supplied along with the read.

decompose :: [Mat44D] -> BamRaw -> Maybe PosPrimChunks
decompose matrices br =
    if isUnmapped b || isDuplicate b || not (isValidRefseq b_rname)
    then Nothing else Just (b_rname, b_pos, pchunks)
  where
    b@BamRec{..} = unpackBam br
    pchunks = firstBase b_pos 0 0 (maybe [] id $ getMd b) matrices

    !max_cig = V.length b_cigar
    !max_seq = V.length b_seq
    !baq     = extAsString "BQ" b

    -- This will compute the effective quality.  As far as I can see
    -- from the BAM spec V1.4, the qualities that matter are QUAL, MAPQ,
    -- and BAQ.  If QUAL is invalid, we replace it (arbitrarily) with
    -- 23 (assuming a rather conservative error rate of ~0.5%), BAQ is
    -- added to QUAL, and MAPQ is an upper limit for effective quality.
    get_seq :: Int -> Nucleotides -> Mat44D -> DamagedBase
    get_seq i = case b_seq `V.unsafeIndex` i of                                 -- nucleotide
            n | n == nucsA -> DB nucA qe
              | n == nucsC -> DB nucC qe
              | n == nucsG -> DB nucG qe
              | n == nucsT -> DB nucT qe
              | otherwise  -> DB nucA (Q 0)
      where
        !q = case b_qual `V.unsafeIndex` i of Q 0xff -> Q 30 ; x -> x           -- quality; invalid (0xff) becomes 30
        !q' | i >= B.length baq = q                                             -- no BAQ available
            | otherwise = Q (unQ q + (B.index baq i - 64))                      -- else correct for BAQ
        !qe = min q' b_mapq                                                     -- use MAPQ as upper limit

    get_seq' :: Int -> Mat44D -> DamagedBase
    get_seq' i = case b_seq `V.unsafeIndex` i of                                -- nucleotide
            n | n == nucsA -> DB nucA qe nucsA
              | n == nucsC -> DB nucC qe nucsC
              | n == nucsG -> DB nucG qe nucsG
              | n == nucsT -> DB nucT qe nucsT
              | otherwise  -> DB nucA (Q 0) n
      where
        !q = case b_qual `V.unsafeIndex` i of Q 0xff -> Q 30 ; x -> x           -- quality; invalid (0xff) becomes 30
        !q' | i >= B.length baq = q                                             -- no BAQ available
            | otherwise = Q (unQ q + (B.index baq i - 64))                      -- else correct for BAQ
        !qe = min q' b_mapq                                                     -- use MAPQ as upper limit

    -- Look for first base following the read's start or a gap (CIGAR
    -- code N).  Indels are skipped, since these are either bugs in the
    -- aligner or the aligner getting rid of essentially unalignable
    -- bases.
    firstBase :: Int -> Int -> Int -> [MdOp] -> [Mat44D] -> PrimChunks
    firstBase !_   !_  !_    _ [        ] = EndOfRead
    firstBase !pos !is !ic mds mms@(m:ms)
        | is >= max_seq || ic >= max_cig = EndOfRead
        | otherwise = case b_cigar `V.unsafeIndex` ic of
            Ins :* cl ->            firstBase  pos (cl+is) (ic+1) mds mms
            SMa :* cl ->            firstBase  pos (cl+is) (ic+1) mds mms
            Del :* cl ->            firstBase (pos+cl) is  (ic+1) (drop_del cl mds) mms
            Nop :* cl ->            firstBase (pos+cl) is  (ic+1) mds mms
            HMa :*  _ ->            firstBase  pos     is  (ic+1) mds mms
            Pad :*  _ ->            firstBase  pos     is  (ic+1) mds mms
            Mat :*  0 ->            firstBase  pos     is  (ic+1) mds mms
            Mat :*  _ -> Seek pos $ nextBase 0 pos     is   ic 0  mds m ms
      where
        -- We have to treat (MdNum 0), because samtools actually
        -- generates(!) it all over the place and if not handled as a
        -- special case, it looks like an incinsistend MD field.
        drop_del n (MdDel ns : mds')
            | n < length ns = MdDel (drop n ns) : mds'
            | n > length ns = drop_del (n - length ns) mds'
            | otherwise     = mds'
        drop_del n (MdNum 0 : mds') = drop_del n mds'
        drop_del _ mds'     = mds'

    -- Generate likelihoods for the next base.  When this gets called,
    -- we are looking at an M CIGAR operation and all the subindices are
    -- valid.
    -- I don't think we can ever get (MdDel []), but then again, who
    -- knows what crazy shit samtools decides to generate.  There is
    -- little harm in special-casing it.
    nextBase :: Int -> Int -> Int -> Int -> Int -> [MdOp] -> Mat44D -> [Mat44D] -> PrimBase
    nextBase !wt !pos !is !ic !io mds m ms = case mds of
        MdNum   0 : mds' -> nextBase wt pos is ic io mds' m ms
        MdDel  [] : mds' -> nextBase wt pos is ic io mds' m ms
        MdNum   1 : mds' -> nextBase' (get_seq' is       m) mds'
        MdNum   n : mds' -> nextBase' (get_seq' is       m) (MdNum (n-1) : mds')
        MdRep ref : mds' -> nextBase' (get_seq  is ref   m) mds'
        MdDel   _ : _    -> nextBase' (get_seq  is nucsN m) mds
        [              ] -> nextBase' (get_seq  is nucsN m) [ ]
      where
        nextBase' ref mds' = Base wt ref b_mapq (isReversed b)
                           $ nextIndel  [] [] (pos+1) (is+1) ic (io+1) mds' ms

    -- Look for the next indel after a base.  We collect all indels (I
    -- and D codes) into one combined operation.  If we hit N or the
    -- read's end, we drop all of it (indels next to a gap indicate
    -- trouble).  Other stuff is skipped: we could check for stuff that
    -- isn't valid in the middle of a read (H and S), but then what
    -- would we do about it anyway?  Just ignoring it is much easier and
    -- arguably at least as correct.
    nextIndel :: [[DamagedBase]] -> [Nucleotides] -> Int -> Int -> Int -> Int -> [MdOp] -> [Mat44D] -> PrimChunks
    nextIndel _   _   !_   !_  !_  !_   _  [        ] = EndOfRead
    nextIndel ins del !pos !is !ic !io mds mms@(m:ms)
        | is >= max_seq || ic >= max_cig = EndOfRead
        | otherwise = case b_cigar `V.unsafeIndex` ic of
            Ins :* cl ->             nextIndel (isq cl) del   pos (cl+is) (ic+1) 0 mds (drop cl mms)
            SMa :* cl ->             nextIndel  ins     del   pos (cl+is) (ic+1) 0 mds (drop cl mms)
            Del :* cl ->             nextIndel ins (del++dsq) (pos+cl) is (ic+1) 0 mds' mms
                where (dsq,mds') = split_del cl mds
            Pad :*  _ ->             nextIndel  ins     del   pos     is  (ic+1) 0 mds mms
            HMa :*  _ ->             nextIndel  ins     del   pos     is  (ic+1) 0 mds mms
            Nop :* cl ->             firstBase               (pos+cl) is  (ic+1)   mds mms  -- ends up generating a 'Seek'
            Mat :* cl | io == cl  -> nextIndel  ins     del   pos     is  (ic+1) 0 mds mms
                      | otherwise -> indel del out $ nextBase (length del) pos is ic io mds m ms -- ends up generating a 'Base'
      where
        indel d o k = rlist o `seq` Indel d o k
        out    = concat $ reverse ins
        isq cl = zipWith ($) [ get_seq i gap | i <- [is..is+cl-1] ] (take cl mms) : ins
        rlist [] = ()
        rlist (a:as) = a `seq` rlist as

        -- We have to treat (MdNum 0), because samtools actually
        -- generates(!) it all over the place and if not handled as a
        -- special case, it looks like an incinsistend MD field.
        split_del n (MdDel ns : mds')
            | n < length ns = (take n ns, MdDel (drop n ns) : mds')
            | n > length ns = let (ns', mds'') = split_del (n - length ns) mds' in (ns++ns', mds'')
            | otherwise     = (ns, mds')
        split_del n (MdNum 0 : mds') = split_del n mds'
        split_del n mds'    = (replicate n nucsN, mds')

-- | Statistics about a genotype call.  Probably only useful for
-- fitlering (so not very useful), but we keep them because it's easy to
-- track them.

data CallStats = CallStats { read_depth       :: {-# UNPACK #-} !Int       -- number of contributing reads
                           , reads_mapq0      :: {-# UNPACK #-} !Int       -- number of (non-)contributing reads with MAPQ==0
                           , sum_mapq         :: {-# UNPACK #-} !Int       -- sum of map qualities of contributing reads
                           , sum_mapq_squared :: {-# UNPACK #-} !Int }     -- sum of squared map qualities of contributing reads
  deriving (Show, Eq)

instance Monoid CallStats where
    mempty      = CallStats { read_depth       = 0
                            , reads_mapq0      = 0
                            , sum_mapq         = 0
                            , sum_mapq_squared = 0 }
    mappend x y = CallStats { read_depth       = read_depth x + read_depth y
                            , reads_mapq0      = reads_mapq0 x + reads_mapq0 y
                            , sum_mapq         = sum_mapq x + sum_mapq y
                            , sum_mapq_squared = sum_mapq_squared x + sum_mapq_squared y }

-- | Genotype likelihood values.  A variant call consists of a position,
-- some measure of qualities, genotype likelihood values, and a
-- representation of variants.  A note about the GL values:  @VCF@ would
-- normalize them so that the smallest one becomes zero.  We do not do
-- that here, since we might want to compare raw values for a model
-- test.  We also store them in a 'Double' to make arithmetics easier.
-- Normalization is appropriate when converting to @VCF@.
--
-- If GL is given, we follow the same order used in VCF:
-- \"the ordering of genotypes for the likelihoods is given by:
-- F(j/k) = (k*(k+1)/2)+j.  In other words, for biallelic sites the
-- ordering is: AA,AB,BB; for triallelic sites the ordering is:
-- AA,AB,BB,AC,BC,CC, etc.\"

type GL = U.Vector Prob

newtype V_Nuc  = V_Nuc  (U.Vector Nucleotide)  deriving (Eq, Ord, Show)
newtype V_Nucs = V_Nucs (U.Vector Nucleotides) deriving (Eq, Ord, Show)

data IndelVariant = IndelVariant { deleted_bases  :: !V_Nucs, inserted_bases :: !V_Nuc }
      deriving (Eq, Ord, Show)


-- | Map quality and a list of encountered bases, with damage
-- information and reference base if known.
type BasePile  = [( Qual,                  DamagedBase   )]

-- | Map quality and a list of encountered indel variants.  The deletion
-- has the reference sequence, if known, an insertion has the inserted
-- sequence with damage information.
type IndelPile = [( Qual, ([Nucleotides], [DamagedBase]) )]   -- a list of indel variants

-- | Running pileup results in a series of piles.  A 'Pile' has the
-- basic statistics of a 'VarCall', but no GL values and a pristine list
-- of variants instead of a proper call.  We emit one pile with two
-- 'BasePile's (one for each strand) and one 'IndelPile' (the one
-- immediately following) at a time.

data Pile' a b = Pile { p_refseq     :: {-# UNPACK #-} !Refseq
                      , p_pos        :: {-# UNPACK #-} !Int
                      , p_snp_stat   :: {-# UNPACK #-} !CallStats
                      , p_snp_pile   :: a
                      , p_indel_stat :: {-# UNPACK #-} !CallStats
                      , p_indel_pile :: b }
  deriving Show

type Pile  = Pile' (BasePile, BasePile) IndelPile

-- | The pileup enumeratee takes 'BamRaw's, decomposes them, interleaves
-- the pieces appropriately, and generates 'Pile's.  The output will
-- contain at most one 'BasePile' and one 'IndelPile' for each position,
-- piles are sorted by position.
--
-- This top level driver receives 'BamRaw's.  Unaligned reads and
-- duplicates are skipped (but not those merely failing quality checks).
-- Processing stops when the first read with invalid 'br_rname' is
-- encountered or a t end of file.

{-# INLINE pileup #-}
pileup :: Enumeratee [PosPrimChunks] [Pile] IO a
pileup = eneeCheckIfDonePass (icont . runPileM pileup' finish (Refseq 0) 0 [] Empty)
  where
    finish () _r _p [] Empty out inp = idone (liftI out) inp
    finish () _ _ _ _ _ _ = error "logic error: leftovers after pileup"


-- | The pileup logic keeps a current coordinate (just two integers) and
-- two running queues: one of /active/ 'PrimBase's that contribute to
-- current genotype calling and on of /waiting/ 'PrimBase's that will
-- contribute at a later point.
--
-- Oppan continuation passing style!  Not only is the CPS version of the
-- state monad (we have five distinct pieces of state) somewhat faster,
-- we also need CPS to interact with the mechanisms of 'Iteratee'.  It
-- makes implementing 'yield', 'peek', and 'bump' straight forward.

newtype PileM m a = PileM { runPileM :: forall r . (a -> PileF m r) -> PileF m r }

-- | The things we drag along in 'PileM'.  Notes:
-- * The /active/ queue is a simple stack.  We add at the front when we
--   encounter reads, which reverses them.  When traversing it, we traverse
--   reads backwards, but since we accumulate the 'BasePile', it gets reversed
--   back.  The new /active/ queue, however, is no longer reversed (as it should
--   be).  So after the traversal, we reverse it again.  (Yes, it is harder to
--   understand than using a proper deque type, but it is cheaper.
--   There may not be much point in the reversing, though.)

type PileF m r = Refseq -> Int ->                               -- current position
                 [PrimBase] ->                                  -- active queue
                 Heap ->                                        -- waiting queue
                 (Stream [Pile] -> Iteratee [Pile] m r) ->      -- output function
                 Stream [PosPrimChunks] ->                      -- pending input
                 Iteratee [PosPrimChunks] m (Iteratee [Pile] m r)

instance Functor (PileM m) where
    {-# INLINE fmap #-}
    fmap f (PileM m) = PileM $ \k -> m (k . f)

instance Applicative (PileM m) where
    {-# INLINE pure #-}
    pure a = PileM $ \k -> k a
    {-# INLINE (<*>) #-}
    u <*> v = PileM $ \k -> runPileM u (\a -> runPileM v (k . a))

instance Monad (PileM m) where
    {-# INLINE return #-}
    return a = PileM $ \k -> k a
    {-# INLINE (>>=) #-}
    m >>=  k = PileM $ \k' -> runPileM m (\a -> runPileM (k a) k')

{-# INLINE get_refseq #-}
get_refseq :: PileM m Refseq
get_refseq = PileM $ \k r -> k r r

{-# INLINE get_pos #-}
get_pos :: PileM m Int
get_pos = PileM $ \k r p -> k p r p

{-# INLINE upd_pos #-}
upd_pos :: (Int -> Int) -> PileM m ()
upd_pos f = PileM $ \k r p -> k () r $! f p

-- | Sends one piece of output downstream.  You are not expected to
-- understand how this works, but inlining 'eneeCheckIfDone' plugged an
-- annoying memory leak.
{-# INLINE yieldPile #-}
yieldPile :: CallStats -> BasePile -> BasePile -> CallStats -> IndelPile -> PileM m ()
yieldPile x1 x2a x2b x3 x4 = PileM $ \ !kont !r !p !a !w !out !inp -> Iteratee $ \od oc ->
      let recurse           = kont () r p a w
          onDone y s        = od (idone y s) inp
          onCont k Nothing  = runIter (recurse k inp) od oc
          onCont k (Just e) = runIter (throwRecoverableErr e (recurse k . (<>) inp)) od oc
          pile              = Pile r p x1 (x2a,x2b) x3 x4
      in runIter (out (Chunk [pile])) onDone onCont

-- | The actual pileup algorithm.  If /active/ contains something,
-- continue here.  Else find the coordinate to continue from, which is
-- the minimum of the next /waiting/ coordinate and the next coordinate
-- in input; if found, continue there, else we're all done.
pileup' :: PileM m ()
pileup' = PileM $ \ !k !refseq !pos !active !waiting !out !inp ->

    let recurse     = runPileM pileup'  k refseq pos active waiting out
        cont2 rs po = runPileM pileup'' k     rs po  active waiting out inp
        leave       =                k () refseq pos active waiting out inp

    in case (active, getMinKeyH waiting, inp) of
        ( _:_,       _,                _  ) -> cont2 refseq pos
        ( [ ], Just nw, EOF            _  ) -> cont2 refseq nw
        ( [ ], Nothing, EOF            _  ) -> leave
        (   _,       _, Chunk [         ] ) -> liftI recurse
        ( [ ], Nothing, Chunk ((r,p,_):_) ) -> cont2 r p
        ( [ ], Just nw, Chunk ((r,p,_):_) )
                     | (refseq,nw) <= (r,p) -> cont2 refseq nw
                     | otherwise            -> cont2 r p

pileup'' :: PileM m ()
pileup'' = do
    -- Input is still 'BamRaw', since these can be relied on to be
    -- sorted.  First see if there is any input at the current location,
    -- if so, decompose it and add it to the appropriate queue.
    p'feed_input
    p'check_waiting
    ((fin_bsL, fin_bpL), (fin_bsR, fin_bpR), (fin_is, fin_ip)) <- p'scan_active

    -- Output, but don't bother emitting empty piles.  Note that a plain
    -- basecall still yields an entry in the 'IndelPile'.  This is necessary,
    -- because actual indel calling will want to know how many reads /did not/
    -- show the variant.  However, if no reads show any variant, and here is the
    -- first place where we notice that, the pile is useless.
    let uninteresting (_,(d,i)) = null d && null i
    unless (null fin_bpL && null fin_bpR && all uninteresting fin_ip) $
        yieldPile (fin_bsL <> fin_bsR) fin_bpL fin_bpR fin_is fin_ip

    -- Bump coordinate and loop.  (Note that the bump to the next
    -- reference /sequence/ is done implicitly, because we will run out of
    -- reads and restart in 'pileup''.)
    upd_pos succ
    pileup'

-- | Feeds input as long as it starts at the current position
p'feed_input :: PileM m ()
p'feed_input = PileM $ \kont rs po ac wt out inp -> case inp of
        Chunk [   ] -> liftI $ runPileM p'feed_input kont rs po ac wt out
        Chunk ((rs', po', prim):bs)
            | rs == rs' && po == po' ->
                case prim of
                    Seek   !p !pb -> let wt' = Node p pb Empty Empty `unionH` wt
                                     in runPileM p'feed_input kont rs po    ac   wt' out (Chunk bs)
                    Indel _ _ !pb ->    runPileM p'feed_input kont rs po (pb:ac) wt  out (Chunk bs)
                    EndOfRead     ->    runPileM p'feed_input kont rs po     ac  wt  out (Chunk bs)
        _           -> kont () rs po ac wt out inp

-- | Checks /waiting/ queue.  If there is anything waiting for the
-- current position, moves it to /active/ queue.
p'check_waiting :: PileM m ()
p'check_waiting = PileM $ \kont rs po ->
        let go ac !wt = case viewMinH wt of
                Just (!mk, !pb, !wt') | mk == po -> go (pb:ac) wt'
                _                                -> kont () rs po ac wt
        in go


-- | Scans /active/ queue and makes a 'BasePile'.  Also sees what's next
-- in the 'PrimChunks':  'Indel's contribute to an 'IndelPile', 'Seek's
-- and deletions are pushed back to the /waiting/ queue, 'EndOfRead's
-- are removed, and everything else is added to the fresh /active/
-- queue.
p'scan_active :: PileM m (( CallStats, [( Qual, DamagedBase )] ),
                          ( CallStats, [( Qual, DamagedBase )] ),
                          ( CallStats, [( Qual, ([Nucleotides], [DamagedBase]) )] ))
p'scan_active =
    PileM $ \kont rs pos ac0 wt -> go (\r -> kont r rs pos) [] wt mempty mempty mempty ac0
  where
    go k !ac !wt !bpileL !bpileR !ipile [                               ] = k (bpileL, bpileR, ipile) (reverse ac) wt
    go k !ac !wt !bpileL !bpileR !ipile (Base nwt qs mq str pchunks : bs) =
        case pchunks of
            _ | nwt > 0 ->        b' `seq` go k  (b':ac)            wt     bpileL    bpileR     ipile  bs
            Seek p' pb' | str           -> go k      ac (ins p' pb' wt) (z bpileL)   bpileR     ipile  bs
                        | otherwise     -> go k      ac (ins p' pb' wt)    bpileL (z bpileR)    ipile  bs
            Indel nd ni pb' | str       -> go k (pb':ac)            wt  (z bpileL)   bpileR  (y ipile) bs
                            | otherwise -> go k (pb':ac)            wt     bpileL (z bpileR) (y ipile) bs
                where y  = put mq (nd,ni)
            EndOfRead       | str       -> go k      ac             wt  (z bpileL)   bpileR     ipile  bs
                            | otherwise -> go k      ac             wt     bpileL (z bpileR)    ipile  bs
        where
            b' = Base (nwt-1) qs mq str pchunks
            z  = put mq qs

    ins q v w = Node q v Empty Empty `unionH` w

    put (Q !q) !x (!st,!vs) = ( st { read_depth       = read_depth st + 1
                                   , reads_mapq0      = reads_mapq0 st + (if q == 0 then 1 else 0)
                                   , sum_mapq         = sum_mapq st + fromIntegral q
                                   , sum_mapq_squared = sum_mapq_squared st + fromIntegral q * fromIntegral q }
                              , (Q q, x) : vs )


-- | We need a simple priority queue.  Here's a skew heap (specialized
-- to strict 'Int' priorities and 'PrimBase' values).
data Heap = Empty | Node {-# UNPACK #-} !Int PrimBase Heap Heap

unionH :: Heap -> Heap -> Heap
Empty                 `unionH` t2                    = t2
t1                    `unionH` Empty                 = t1
t1@(Node k1 x1 l1 r1) `unionH` t2@(Node k2 x2 l2 r2)
   | k1 <= k2                                       = Node k1 x1 (t2 `unionH` r1) l1
   | otherwise                                      = Node k2 x2 (t1 `unionH` r2) l2

getMinKeyH :: Heap -> Maybe Int
getMinKeyH Empty          = Nothing
getMinKeyH (Node x _ _ _) = Just x

viewMinH :: Heap -> Maybe (Int, PrimBase, Heap)
viewMinH Empty          = Nothing
viewMinH (Node k v l r) = Just (k, v, l `unionH` r)

