{-# LANGUAGE BangPatterns, Rank2Types #-}
{-# OPTIONS_GHC -funbox-strict-fields #-}
module Bio.Bam.Pileup where

import Debug.Trace
import Text.Printf

import Bio.Base
import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Iteratee
import Bio.Util ( (<#>) )

import Control.Applicative
import Control.Monad hiding ( mapM_ )
import Control.Monad.Fix ( fix )
import Data.Foldable
import Data.List ( tails, intercalate, intersperse )
import Numeric ( showFFloat )

import qualified Data.ByteString        as B
import qualified Data.Vector.Unboxed    as V
import qualified Data.Set               as Set

import Prelude hiding ( foldr, concat, mapM_, all, sum )

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
-- For aDNA, we need a substitution probability.  We have three options:
-- use an empirically determined PSSM, use an arithmetically defined
-- PSSM based on the /Johnson/ model, use a context sensitive PSSM based
-- on the /Johnson/ model and an alignment.  Using /Dindel/, actual
-- substitutions relative to a called haplotype would be taken into
-- account.  Since we're not going to do that, taking alignments into
-- account is difficult, somewhat approximate, and therefore not worth
-- the hassle.
--
-- Regarding the error model, there's a choice between /samtools/maq/ or
-- the naive model everybody else uses.  Naive is easy to marry to aDNA,
-- samtools is (probably) better.  Either way, we introduce a number of
-- parameters (@eta@ and @kappa@ for /samtools/, @lambda@, @delta@,
-- @delta_ss@ for /Johnson/).  Running a maximum likehood fit for those
-- may be valuable.  It would be cool, if we could do that without
-- rerunning the complete genotype caller, but it's not a priority.
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


-- *TODO*
--
-- * A whole lot of testing.
-- * Actual genotype calling.
-- * ML fitting and evaluation of parameters for different possible
--   error and damage models.
-- * Check the 'decompose' logic, in particular, make sure the waiting
--   time after deletions is exactly right
-- * Ploidy must be a variable; we definitely need calling for diploid
--   and haploid genomes; higher ploidy is nice to have if it comes at
--   acceptable cost.

-- | For likelihoods for bases @A, C, G, T@.
data Double4 = D4 !Double !Double !Double !Double

instance Show Double4 where
    showsPrec _ (D4 a c g t) = foldr (.) id . intersperse ((:) ' ') $ map (showFFloat (Just 2)) [a,c,g,t]

-- | The primitive pieces for genotype calling:  A position, a base
-- represented as four likelihoods, an inserted sequence, and the
-- length of a deleted sequence.  The logic is that we look at a base
-- followed by some indel, and all those indels are combined into a
-- single insertion and a single deletion.
data PrimChunks = Seek !Int PrimBase                            -- ^ skip to position (at start or after N operation)
                | Indel !Int [(Nucleotide, Qual)] PrimBase      -- ^ observed deletion and insertion between two bases
                | EndOfRead                                     -- ^ nothing anymore
  deriving Show

data PrimBase = Base { _pb_wait   :: !Int                       -- ^ number of bases to wait due to a deletion
                     , _pb_likes  :: !Double4                   -- ^ four likelihoods
                     , _pb_mapq   :: !Qual                      -- ^ map quality
                     , _pb_chunks :: PrimChunks }               -- ^ more chunks
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
        !q' | i >= B.length baq = q                                 -- no BAQ available
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
            (Mat, _) -> Seek pos $ nextBase 0 pos     is   ic 0


    -- Generate likelihoods for the next base.  When this gets called,
    -- we are looking at an M CIGAR operation and all the subindices are
    -- valid.
    nextBase :: Int -> Int -> Int -> Int -> Int -> PrimBase
    nextBase !wt !pos !is !ic !io = Base wt (get_seq (dm is) is) mapq $ nextIndel  [] 0 (pos+1) (is+1) ic (io+1)


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
                     | otherwise -> Indel del out $ nextBase del pos is   ic   io   -- ends up generating a 'Base'
            (Nop,cl) ->             firstBase               (pos+cl) is  (ic+1)     -- ends up generating a 'Seek'
      where
        out    = concat $ reverse ins
        isq cl = [ get_seq (,) i | i <- [is..is+cl-1] ] : ins


-- | A 'DamageModel' is a function that gives likelihoods for all
-- possible four bases, given the sequenced base, the quality, and the
-- position in the read.  That means it can't take the actual alignment
-- into account... but nobody seems too keen on doing that anyway.
type DamageModel = Int              -- ^ position in read
                -> Nucleotide       -- ^ base
                -> Qual             -- ^ quality score
                -> Double4          -- ^ results in four likelihoods

-- | 'DamageModel' for undamaged DNA.  The likelihoods follow directly
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
-- genotype likelihood values, and a representation of variants.  A note
-- about the 'vc_pl' values:  @VCF@ would normalize them so that the
-- smallest one becomes zero.  We do not do that here, since we might
-- want to compare raw values for a model test.  We also store them in a
-- 'Double' to make arithmetics easier.  Normalization is appropriate
-- when converting to @VCF@.  Also note that 'vc_pl' can be empty at the
-- stage where we collected variants, but did not do proper var calling.
--
-- If 'vc_pl' is given, we follow the same order used in VCF:
-- "the ordering of genotypes for the likelihoods is given by:
-- F(j/k) = (k*(k+1)/2)+j.  In other words, for biallelic sites the
-- ordering is: AA,AB,BB; for triallelic sites the ordering is:
-- AA,AB,BB,AC,BC,CC, etc."

data VarCall a = VarCall { vc_refseq     :: !Refseq
                         , vc_pos        :: !Int
                         , vc_depth      :: !Int                  -- number of contributing reads
                         , vc_mapq0      :: !Int                  -- number of contributing reads with MAPQ==0
                         , vc_sum_mapq   :: !Int                  -- sum of map qualities of contributring reads
                         , vc_sum_mapq2  :: !Int                  -- sum of squared map qualities of contributing reads
                         , vc_pl         :: !(V.Vector Double)    -- PL values in dB
                         , vc_vars       :: a }                   -- variant calls, depending on context

type SnpVars = ()                           -- no additonal info needed for SNP calls
type IndelVars = [( Int, [Nucleotide] )]    -- indel variant: number of deletions, inserted sequence

-- Both types of piles carry along the map quality.  We'll only need it
-- in the case of Indels.
type BasePile  = [( Qual, Double4 )]                       -- a list of encountered bases
type IndelPile = [( Qual, (Int, [(Nucleotide, Qual)]) )]   -- a list of indel variants

-- | Running pileup results in a series of piles.  A 'Pile' has the
-- basic statistics of a 'VarCall', but no PL values and a pristine list
-- of variants instead of a proper call.

type Pile = Either (VarCall BasePile) (VarCall IndelPile)


-- | The pileup enumeratee takes 'BamRaw's, decomposes them, interleaves
-- the pieces appropriately, and generates 'Pile's.  The output will
-- contain at most one 'BasePile' and one 'IndelPile' for each position,
-- piles are sorted by position.
--
-- This top level driver receives 'BamRaw's.  Unaligned reads and
-- duplicates are skipped (but not those merely failing quality checks).
-- Processing stops when the first read with invalid 'br_rname' is
-- encountered or a t end of file.

pileup :: MonadIO m => DamageModel -> Enumeratee [BamRaw] [Pile] m a
pileup dm = takeWhileE (isValidRefseq . br_rname) ><> filterStream useable ><>
            eneeCheckIfDone (liftI . runPileM pileup' finish (Refseq 0) 0 [] Empty dm)
  where
    useable br = not (br_isUnmapped br || br_isDuplicate br)

    finish () _r _p [] Empty _dm out inp = idone (liftI out) inp
    finish () _ _ _ _ _ _ _ = error "logic error: leftovers after pileup"


-- | The pileup logic keeps a current coordinate (just two integers) and
-- two running queues: one of /active/ 'PrimBases' that contribute to
-- current genotype calling and on of /waiting/ 'PrimBases' that will
-- contribute at later point.
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
--   understand than using a proper deque type, but it is cheaper.)

type PileF m r = Refseq -> Int ->                               -- current position
                 [PrimBase] ->                                  -- active queue
                 Heap PrimBase ->                               -- waiting queue
                 DamageModel ->
                 (Stream [Pile] -> Iteratee [Pile] m r) ->      -- output function
                 Stream [BamRaw] ->                             -- pending input
                 Iteratee [BamRaw] m (Iteratee [Pile] m r)

instance Functor (PileM m) where
    fmap f (PileM m) = PileM $ \k -> m (k . f)

instance Monad (PileM m) where
    return a = PileM $ \k -> k a
    m >>=  k = PileM $ \k' -> runPileM m (\a -> runPileM (k a) k')

instance MonadIO m => MonadIO (PileM m) where
    liftIO m = PileM $ \k r p a w d o i -> liftIO m >>= \x -> k x r p a w d o i

get_refseq :: PileM m Refseq
get_refseq = PileM $ \k r -> k r r

get_pos :: PileM m Int
get_pos = PileM $ \k r p -> k p r p

upd_pos :: (Int -> Int) -> PileM m ()
upd_pos f = PileM $ \k r p -> k () r $! f p

set_pos :: (Refseq, Int) -> PileM m ()
set_pos (!r,!p) = PileM $ \k _ _ -> k () r p

get_active :: PileM m [PrimBase]
get_active = PileM $ \k r p a -> k a r p a

upd_active :: ([PrimBase] -> [PrimBase]) -> PileM m ()
upd_active f = PileM $ \k r p a -> k () r p $! f a

get_waiting :: PileM m (Heap PrimBase)
get_waiting = PileM $ \k r p a w -> k w r p a w

upd_waiting :: (Heap PrimBase -> Heap PrimBase) -> PileM m ()
upd_waiting f = PileM $ \k r p a w -> k () r p a $! f w

get_damage_model :: PileM m DamageModel
get_damage_model = PileM $ \k r p a w d -> k d r p a w d

maybe_min :: Ord a => Maybe a -> Maybe a -> Maybe a
maybe_min Nothing   Nothing = Nothing
maybe_min (Just a)  Nothing = Just a
maybe_min Nothing  (Just b) = Just b
maybe_min (Just a) (Just b) = Just $ min a b

yield :: Monad m => Pile -> PileM m ()
yield x = PileM $ \k r p a w d out inp ->
    eneeCheckIfDone (\out' -> k () r p a w d out' inp) . out $ Chunk [x]

-- | Inspect next input element, if any.  Returns @Just b@ if @b@ is the
-- next input element, @Nothing@ if no such element exists.  Waits for
-- more input if nothing is available immediately.
peek :: PileM m (Maybe BamRaw)
peek = PileM $ \k r p a w d out inp -> case inp of
        EOF     _   -> k Nothing r p a w d out inp
        Chunk [   ] -> liftI $ runPileM peek k r p a w d out
        Chunk (b:_) -> k (Just b) r p a w d out inp

-- | Discard next input element, if any.  Does nothing if input has
-- already ended.  Waits for input to discard if nothing is available
-- immediately.
bump :: PileM m ()
bump = PileM $ \k r p a w d out inp -> case inp of
        EOF     _   -> k () r p a w d out inp
        Chunk [   ] -> liftI $ runPileM bump k r p a w d out
        Chunk (_:x) -> k () r p a w d out (Chunk x)


consume_active :: a -> (a -> PrimBase -> PileM m a) -> PileM m a
consume_active nil cons = do ac <- get_active
                             upd_active (const [])
                             foldM cons nil ac

-- | The actual pileup algorithm.
pileup' :: MonadIO m => PileM m ()
pileup' = do
    refseq       <- get_refseq
    active       <- get_active
    next_waiting <- fmap ((,) refseq) . getMinKey <$> get_waiting
    next_input   <- fmap (\b -> (br_rname b, br_pos b)) <$> peek

    -- posn         <- get_pos
    -- liftIO $ printf "pileup' @%d:%d, %d active, %d waiting\n"
        -- (unRefseq refseq) posn (length active) (-1::Int)


    -- If /active/ contains something, continue here.  Else find the coordinate
    -- to continue from, which is the minimum of the next /waiting/ coordinate
    -- and the next coordinate in input; if found, continue there, else we're
    -- all done.
    case (active, maybe_min next_waiting next_input) of
        ( (_:_),       _ ) -> pileup''
        ( [   ], Just mp ) -> set_pos mp >> pileup''
        ( [   ], Nothing ) -> return ()

pileup'' :: MonadIO m => PileM m ()
pileup'' = do
    -- Input is still 'BamRaw', since these can be relied on to be
    -- sorted.  First see if there is any input at the current location,
    -- if so, decompose it and add it to the appropriate queue.
    rs <- get_refseq
    po <- get_pos
    dm <- get_damage_model

    -- liftIO $ printf "pileup' @%d:%d, %d active, %d waiting\n"
        -- (unRefseq rs) po (-1::Int) (-1::Int)

    fix $ \loop -> peek >>= mapM_ (\br ->
            when (br_rname br == rs && br_pos br == po) $ do
                bump
                case decompose dm br of
                    Seek    p pb -> upd_waiting (insert p pb)
                    Indel _ _ pb -> upd_active (pb:)
                    EndOfRead    -> return ()
                loop)


    -- Check /waiting/ queue.  If there is anything waiting for the
    -- current position, move it to /active/ queue.
    fix $ \loop -> (getMin <$> get_waiting) >>= mapM_ (\(mk,pb) ->
            when (mk == po) $ do
                upd_active (pb:)
                upd_waiting $ dropMin
                loop)

    -- Scan /active/ queue and make a 'BasePile'.  Also see what's next in the
    -- 'PrimChunks':  'Indel's contribute to an 'IndelPile', 'Seek's and
    -- deletions are pushed back to the /waiting/ queue, 'EndOfRead's are
    -- removed, and everything else is added to the fresh /active/ queue.
    let pile0 = VarCall rs po 0 0 0 0 V.empty
    (fin_bp, fin_ip) <- consume_active (pile0 [], pile0 []) $
        \(!bpile, !ipile) (Base wt qs mq pchunks) ->
                let put (Q q) x vc = vc { vc_depth     = vc_depth vc + 1
                                        , vc_mapq0     = vc_mapq0 vc + (if q == 0 then 1 else 0)
                                        , vc_sum_mapq  = vc_sum_mapq  vc + fromIntegral q
                                        , vc_sum_mapq2 = vc_sum_mapq2 vc + fromIntegral q * fromIntegral q
                                        , vc_vars      = (Q q, x) : vc_vars vc }
                    b' = Base (wt-1) qs mq pchunks
                in case pchunks of
                    _ | wt > 0        -> do upd_active  (b'  :)         ; return (           bpile,                  ipile )
                    Seek p' pb'       -> do upd_waiting (insert p' pb') ; return ( put mq qs bpile,                  ipile )
                    Indel del ins pb' -> do upd_active  (pb' :)         ; return ( put mq qs bpile, put mq (del,ins) ipile )
                    EndOfRead         -> do                               return ( put mq qs bpile,                  ipile )

    -- We just reversed /active/ inplicitly, which is no desaster, but may come
    -- as a surprise downstream.  So reverse it back.
    upd_active reverse

    -- Output, but don't bother emitting empty piles.  Note that a plain
    -- basecall still yields an entry in the 'IndelPile'.  This is necessary,
    -- because actual indel calling will want to know how many reads /did not/
    -- show the variant.  However, if no reads show any variant, and here is the
    -- first place where we notice that, the pile is useless.
    let uninteresting (_,(d,i)) = d == 0 && null i
    unless (null              (vc_vars fin_bp)) $ yield $ Left  fin_bp
    unless (all uninteresting (vc_vars fin_ip)) $ yield $ Right fin_ip

    -- Bump coordinate and loop.
    upd_pos succ
    pileup'


-- | We need a simple priority queue.  Here's a skew heap (lightly
-- specialized to strict 'Int' priorities).
data Heap a = Empty | Node !Int a (Heap a) (Heap a)

union :: Heap a -> Heap a -> Heap a
Empty                 `union` t2                    = t2
t1                    `union` Empty                 = t1
t1@(Node k1 x1 l1 r1) `union` t2@(Node k2 x2 l2 r2)
   | k1 <= k2                                       = Node k1 x1 (t2 `union` r1) l1
   | otherwise                                      = Node k2 x2 (t1 `union` r2) l2

insert :: Int -> a -> Heap a -> Heap a
insert k v heap = Node k v Empty Empty `union` heap

getMinKey :: Heap a -> Maybe Int
getMinKey Empty          = Nothing
getMinKey (Node x _ _ _) = Just x

getMinVal :: Heap a -> Maybe a
getMinVal Empty          = Nothing
getMinVal (Node _ x _ _) = Just x

getMin :: Heap a -> Maybe (Int,a)
getMin Empty          = Nothing
getMin (Node k v _ _) = Just (k,v)

dropMin :: Heap a -> Heap a
dropMin Empty          = error "dropMin on empty queue... are you sure?!"
dropMin (Node _ _ l r) = l `union` r


-- ------------------------------------------------------------------------------------------------------------

-- | A gentoype caller receives a list of bases at a column and is
-- supposed to turn them into something... PL values come to mind.  XXX
-- Trouble is, the 'something' may well change with the concrete caller.
-- (This used to have a random component, which made everything
-- complicated.  In hindsight, that's silly.  Application of this
-- nonsense can be supported externally.)

-- type Conscall = [Double4] -> Double10

-- | Simple indel calling.  We don't bother with it too much, so here's
-- the gist:  We collect variants (simply different variants, details
-- don't matter), so @n@ variants give rise to (n+1)*n/2 PL values.
-- (That's two out of @(n+1)@, the reference allele, represented here as
-- no deletion and no insertion, is there, too.)  To assign these, we
-- need a likelihood for an observed variant given an assumed genotype.
--
-- For variants of equal length, the likelihood is the sum of qualities
-- of mismatching bases, but no higher than the mapping quality.  That
-- is roughly the likelihood of getting the observed sequence even
-- though the real sequence is a different variant.  For variants of
-- different length, the likelihood is the map quality.  This
-- corresponds to the assumption that indel errors in sequencing are
-- much less likely than mapping errors.  Since this hardly our
-- priority, the approximations are declared good enough.
--
-- About the actual computation:  each of the PL values is the product
-- of the corresponding PL value for each input read.  The PL value for
-- an input read is the likehood of getting that read from a

simple_indel_call :: Int -> VarCall IndelPile -> VarCall IndelVars
simple_indel_call ploidy ip = ip { vc_pl = simple_call ploidy mkpls (vc_vars ip), vc_vars = vars' }
  where
    vars' = Set.toList . Set.fromList $ [ (d, map fst i) | (_q,(d,i)) <- vc_vars ip ]
    match = zipWith (\(n,qn) nr -> if n == nr then 0 else fromIntegral $ unQ qn)
    mkpls (q,(d,i)) = let !q' = fromIntegral (unQ q)
                      in [ if d /= dr || length i /= length ir
                           then q' else q' <#> sum (match i ir) | (dr,ir) <- vars' ]

-- | Compute @PL@ values for the simple case.  The simple case is where
-- where we sample 'ploidy' alleles with equal probability and assume
-- that error occur independently from each other.
--
-- The argument 'pls' is a function that computes the likelihood for
-- getting the current read, for every variant assuming that variant was
-- sampled.

simple_call :: Int -> (a -> [Double]) -> [a] -> V.Vector Double
simple_call ploidy pls = foldl1' (V.zipWith (+)) . map step
  where
    foldl1' _ [     ] = V.singleton 0
    foldl1' f (!a:as) = foldl' f a as

    -- "For biallelic sites the ordering is: AA,AB,BB; for triallelic
    -- sites the ordering is: AA,AB,BB,AC,BC,CC, etc."
    --
    -- To get the order right, we reverse the list, recurse, and reverse
    -- the result again.
    step = V.fromList . reverse . mk_pls ploidy . reverse . pls

    mk_pls 0  _ = return (1/0)
    mk_pls n ls = do ls'@(hd:_) <- tails ls
                     (<#>) hd <$> mk_pls (n-1) ls'

smoke_test :: IO ()
smoke_test =
    decodeAnyBamFile "/mnt/scratch/udo/test.bam" >=> run $ \_hdr ->
    joinI $ pileup noDamage $ mapStreamM_ call_and_print
  where
    call_and_print (Right ic) = put . showCall show_indels . simple_indel_call 2 $ ic
    call_and_print (Left  bc) = return () -- put . showCall show_bases                        $ bc

    put f = putStr $ f "\n"

    show_bases :: a -> ShowS
    show_bases _ = (++) "<<bases>>"

    show_indels :: IndelVars -> ShowS
    show_indels = (++) . intercalate "," . map show_indel

    show_indel :: (Int, [Nucleotide]) -> String
    show_indel (d, ins) = shows ins $ '-' : show d


showCall :: (a -> ShowS) -> VarCall a -> ShowS
showCall f vc = shows (vc_refseq vc) . (:) ':' .
                shows (vc_pos vc) . (:) '\t' .
                f (vc_vars vc) . (++) "\tDP=" .
                shows (vc_depth vc) . (++) ":MQ0=" .
                shows (vc_mapq0 vc) . (++) ":MAPQ=" .
                shows mapq . (:) '\t' .
                show_pl (vc_pl vc)
  where
    show_pl :: V.Vector Double -> ShowS
    show_pl = (++) . intercalate "," . map show . V.toList

    mapq = vc_sum_mapq vc `div` vc_depth vc
