{-# LANGUAGE BangPatterns, RecordWildCards, FlexibleContexts, DeriveGeneric #-}
module Bio.Adna (
    DmgStats(..),
    CompositionStats,
    SubstitutionStats,
    damagePatternsIter,
    damagePatternsIterMD,
    damagePatternsIter2Bit,

    DamageParameters(..),
    NewDamageParameters(..),
    GenDamageParameters(..),
    DamageModel,
    bang,
    Alignment(..),
    FragType(..),
    Subst(..),
    NPair,

    noDamage,
    univDamage,
    empDamage,
    Mat44D(..),
    scalarMat
                ) where

import Bio.Bam
import Bio.Prelude
import Bio.TwoBit
import Bio.Util.Pretty

import qualified Data.Vector.Generic            as G
import qualified Data.Vector.Storable           as VS
import qualified Data.Vector.Unboxed            as U
import qualified Data.Vector.Unboxed.Mutable    as UM

-- ^ Things specific to ancient DNA, e.g. damage models.
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
-- We represent substitution matrices by the type 'Mat44D'.  Internally,
-- this is a vector of packed vectors.  Conveniently, each of the packed
-- vectors represents all transition /into/ the given nucleotide.

newtype Mat44D = Mat44D (U.Vector Double) deriving Show

-- | A 'DamageModel' is a function that gives substitution matrices for
-- each position in a read.  The 'DamageModel' can depend on whether the
-- alignment is reversed, the length of the read and the position.  (In
-- practice, we should probably memoize precomputed damage models
-- somehow.)

type DamageModel = Bool -> Int -> Int -> Mat44D
data Subst = Nucleotide :-> Nucleotide deriving (Eq, Ord, Ix, Show)

infix 9 :->
infix 8 `bang`

-- | Convenience function to access a substitution matrix that has a
-- mnemonic reading.
{-# INLINE bang #-}
bang :: Mat44D -> Subst -> Double
bang (Mat44D v) (N x :-> N y) = v U.! (fromIntegral x + 4 * fromIntegral y)

scalarMat :: Double -> Mat44D
scalarMat s = Mat44D $ U.fromListN 16 [ s, 0, 0, 0
                                      , 0, s, 0, 0
                                      , 0, 0, s, 0
                                      , 0, 0, 0, s ]

-- | 'DamageModel' for undamaged DNA.  The likelihoods follow directly
-- from the quality score.  This needs elaboration to see what to do
-- with amibiguity codes (even though those haven't actually been
-- observed in the wild).

noDamage :: DamageModel
noDamage _ _ _ = scalarMat 1


-- | Parameters for the universal damage model.
--
-- We assume the correct model is either no damage, or single strand
-- damage, or double strand damage.  Each of them comes with a
-- probability.  It turns out that blending them into one is simply
-- accomplished by multiplying these probabilities onto the deamination
-- probabilities.
--
-- For single stranded library prep, only one kind of damage occurs (C
-- to T), it occurs at low frequency ('ssd_delta') everywhere, at high
-- frequency ('ssd_sigma') in single stranded parts, and the overhang
-- length is distributed exponentially with parameter 'ssd_lambda' at
-- the 5' end and 'ssd_kappa' at the 3' end.  (Without UDG treatment,
-- those will be equal.  With UDG, those are much smaller and in fact
-- don't literally represent overhangs.)
--
-- For double stranded library prep, we get C->T damage at the 5' end
-- and G->A at the 3' end with rate 'dsd_sigma' and both in the interior
-- with rate 'dsd_delta'.  Everything is symmetric, and therefore the
-- orientation of the aligned read doesn't matter either.  Both
-- overhangs follow a distribution with parameter 'dsd_lambda'.

data DamageParameters float = DP { ssd_sigma  :: !float         -- deamination rate in ss DNA, SS model
                                 , ssd_delta  :: !float         -- deamination rate in ds DNA, SS model
                                 , ssd_lambda :: !float         -- param for geom. distribution, 5' end, SS model
                                 , ssd_kappa  :: !float         -- param for geom. distribution, 3' end, SS model
                                 , dsd_sigma  :: !float         -- deamination rate in ss DNA, DS model
                                 , dsd_delta  :: !float         -- deamination rate in ds DNA, DS model
                                 , dsd_lambda :: !float }       -- param for geom. distribution, DS model
  deriving (Read, Show, Generic)

instance Pretty f => Pretty (DamageParameters f)
instance Parse f  => Parse  (DamageParameters f)

data NewDamageParameters vec float = NDP { dp_gc_frac :: !float
                                         , dp_mu      :: !float
                                         , dp_nu      :: !float
                                         , dp_alpha5  :: !(vec float)
                                         , dp_beta5   :: !(vec float)
                                         , dp_alpha   :: !float
                                         , dp_beta    :: !float
                                         , dp_alpha3  :: !(vec float)
                                         , dp_beta3   :: !(vec float) }
  deriving (Read, Show, Generic)

instance (Pretty (v f), Pretty f) => Pretty (NewDamageParameters v f)
instance (Parse  (v f), Parse  f) => Parse  (NewDamageParameters v f)

data GenDamageParameters vec float
    = UnknownDamage
    | OldDamage (DamageParameters float)
    | NewDamage (NewDamageParameters vec float)
  deriving (Show, Generic, Read)

instance (Pretty (v f), Pretty f) => Pretty (GenDamageParameters v f)
instance (Parse  (v f), Parse  f) => Parse  (GenDamageParameters v f)




-- | Generic substitution matrix, has C->T and G->A deamination as
-- parameters.  Setting 'p' or 'q' to 0 as appropriate makes this apply
-- to the single stranded or undamaged case.

{-# INLINE genSubstMat #-}
genSubstMat :: Double -> Double -> Mat44D
genSubstMat p q = Mat44D $ U.fromListN 16 [ 1,  0,   q,  0
                                          , 0, 1-p,  0,  0
                                          , 0,  0,  1-q, 0
                                          , 0,  p,   0,  1 ]

univDamage :: DamageParameters Double -> DamageModel
univDamage DP{..} r l i = genSubstMat (p1+p2) (q1+q2)
    where
        (p1, q1) = if r then let lam5 = ssd_lambda ^ (l-i)
                                 lam3 = ssd_kappa ^ (1+i)
                                 lam  = lam3 + lam5 - lam3 * lam5
                                 p    = ssd_sigma * lam + ssd_delta * (1-lam)
                             in (0,p)
                        else let lam5 = ssd_lambda ^ (1+i)
                                 lam3 = ssd_kappa ^ (l-i)
                                 lam  = lam3 + lam5 - lam3 * lam5
                                 p    = ssd_sigma * lam + ssd_delta * (1-lam)
                             in (p,0)

        p2      = dsd_sigma * lam5_ds + dsd_delta * (1-lam5_ds)
        q2      = dsd_sigma * lam3_ds + dsd_delta * (1-lam3_ds)
        lam5_ds = dsd_lambda ^ (1+i)
        lam3_ds = dsd_lambda ^ (l-i)

{-# SPECIALIZE empDamage :: NewDamageParameters U.Vector Double -> DamageModel #-}
empDamage :: G.Vector v Double => NewDamageParameters v Double -> DamageModel
empDamage NDP{..} r l
    | r         = get (flip genSubstMat')
    | otherwise = get       genSubstMat'
  where
    get k i | i+i  <  l = k (fromMaybe dp_alpha (dp_alpha5 G.!? i))
                            (fromMaybe dp_beta  (dp_beta5  G.!? i))
            | otherwise = k (fromMaybe dp_alpha (dp_alpha3 G.!? (l-i-1)))
                            (fromMaybe dp_beta  (dp_beta3  G.!? (l-i-1)))

    genSubstMat' a b = genSubstMat (recip $ 1 + exp (-a)) (recip $ 1 + exp (-b))

-- | Collected \"traditional\" statistics:
--
-- * Base composition near 5' end and near 3' end.  Each consists of
--   five vectors of counts of A,C,G,T, and everything else.
--   'basecompo5' begins with 'context' bases to the left of the 5' end,
--   'basecompo3' ends with 'context' bases to the right of the 3' end.
--
-- * Substitutions.  Counted from the reconstructed alignment, once
--   around the 5' end and once around the 3' end.  For a total of 2*4*4
--   different substitutions.  Positions where the query has a gap are
--   skipped.
--
-- * Substitutions at CpG motifs.  Also counted from the reconstructed
--   alignment, and a CpG site is simply the sequence CG in the
--   reference.  Gaps may confuse that definition, so that CpHpG still
--   counts as CpG, because the H is gapped.  That might actually
--   be desirable.
--
-- * Conditional substitutions.  The 5' and 3' ends count as damaged if
--   the very last position has a C-to-T substitution.  With that in
--   mind, 'substs5d5', 'substs5d3', 'substs5dd' are like 'substs5', but
--   counting only reads where the 5' end is damaged, where the 3' end
--   is damaged, and where both ends are damaged, respectively.
--
-- XXX  This got kind of ugly.  We'll see where this goes...

data DmgStats a = DmgStats {
    basecompo5 :: CompositionStats,
    basecompo3 :: CompositionStats,
    substs5    :: SubstitutionStats,
    substs3    :: SubstitutionStats,
    substs5d5  :: SubstitutionStats,
    substs3d5  :: SubstitutionStats,
    substs5d3  :: SubstitutionStats,
    substs3d3  :: SubstitutionStats,
    substs5dd  :: SubstitutionStats,
    substs3dd  :: SubstitutionStats,
    substs5cpg :: SubstitutionStats,
    substs3cpg :: SubstitutionStats,
    stats_more :: a }

type CompositionStats  = [( Maybe Nucleotide, U.Vector Int )]
type SubstitutionStats = [( Subst, U.Vector Int )]


data FragType = Complete | Leading | Trailing deriving (Show, Eq)
type NPair = ( Nucleotides, Nucleotides )

-- Alignment record, might have been gotten from practically anywhere
-- with varying completeness.  We record anything we can get, most is
-- optional.  Reference sequence is filled with Ns if missing.
data Alignment = ALN
    { a_sequence :: !(U.Vector NPair)       -- the alignment proper
    , a_fragment_type :: !FragType }    -- was the adapter trimmed?


-- | Enumeratee (almost) that computes some statistics from plain BAM
-- (no MD field needed) and a 2bit file.  The 'Alignment' is also
-- reconstructed and passed downstream.  The result of any downstream
-- processing is available in the 'stats_more' field of the result.
--
-- * Get the reference sequence including both contexts once.  If this
--   includes invalid sequence (negative coordinate), pad suitably.
-- * Accumulate counts for the valid parts around 5' and 3' ends as
--   appropriate from flags and config.
-- * Combine the part that was aligned to (so no context) with the read
--   to reconstruct the alignment.
--
-- Arguments are the table of reference names, the 2bit file with the
-- reference, the amount of context outside the alignment desired, the
-- amount of context inside desired, and whether to compensate for
-- leeHom breakage.
--
-- For 'Complete' fragments, we cut the read in the middle, so the 5'
-- and 3' plots stay clean from each other's influence.  'Leading' and
-- 'Trailing' fragments count completely towards the appropriate end.

damagePatternsIter2Bit :: MonadIO m
                       => Refs -> TwoBitFile -> Int -> Int -> Bool
                       -> Iteratee [Alignment] m b
                       -> Iteratee [BamRaw] m (DmgStats b)
damagePatternsIter2Bit refs tbf ctx rng leeHom it =
    mapMaybeStream (\br -> do
        let b@BamRec{..} = unpackBam br
        guard (not $ isUnmapped b)
        let ref_nm = sq_name $ getRef refs b_rname
            ref    = getFragment tbf ref_nm (b_pos - ctx) (alignedLength b_cigar + 2*ctx)
            pairs  = aln_from_ref (U.drop ctx ref) b_seq b_cigar
        return (b, ref, pairs)) =$
    damagePatternsIter ctx rng leeHom it

-- | Enumeratee (almost) that computes some statistics from plain BAM
-- with a valid MD field.  The 'Alignment' is also reconstructed and
-- passed downstream.  The result of any downstream processing is
-- available in the 'stats_more' field of the result.
--
-- * Reconstruct the alignment from CIGAR, SEQ, and MD.
-- * Filter the alignment to get the reference sequence, accumulate it.
-- * Accumulate everything over the alignment.
--
-- Arguments are just the amount of context inside desired.
-- Arguments are the amount of context inside desired, and whether to
-- compensate for leeHom breakage.
--
-- For 'Complete' fragments, we cut the read in the middle, so the 5'
-- and 3' plots stay clean from each other's influence.  'Leading' and
-- 'Trailing' fragments count completely towards the appropriate end.

damagePatternsIterMD :: MonadIO m
                     => Int -> Bool
                     -> Iteratee [Alignment] m b
                     -> Iteratee [BamRaw] m (DmgStats b)
damagePatternsIterMD rng leeHom it =
    mapMaybeStream (\br -> do
        let b@BamRec{..} = unpackBam br
        guard (not $ isUnmapped b)
        md <- getMd b
        let pairs = aln_from_md b_seq b_cigar md
            ref   = U.map fromN $ U.filter ((/=) gap . fst) pairs
        return (b, ref, pairs)) =$
    damagePatternsIter 0 rng leeHom it
  where
    fromN (ns,_) | ns == nucsA = 2
                 | ns == nucsC = 1
                 | ns == nucsG = 3
                 | ns == nucsT = 0
                 | otherwise   = 4

-- | Common logic for statistics. The function 'get_ref_and_aln'
-- reconstructs reference sequence and alignment from a Bam record.  It
-- is expected to construct the alignment with respect to the forwards
-- strand of the reference; we reverse-complement it if necessary.
damagePatternsIter :: MonadIO m
                   => Int -> Int
                   -> Bool -> Iteratee [Alignment] m b
                   -> Iteratee [(BamRec, U.Vector Word8, U.Vector NPair)] m (DmgStats b)
damagePatternsIter ctx rng leeHom it = mapStream revcom_both =$ do
    let maxwidth = ctx + rng
    acc_bc <- liftIO $ UM.replicate (2 * 5 *    maxwidth) (0::Int)
    acc_st <- liftIO $ UM.replicate (2 * 4 * 4 * 4 * rng) (0::Int)
    acc_cg <- liftIO $ UM.replicate (2 * 2 * 4 *     rng) (0::Int)

    it' <- flip mapStreamM it $ \(b@BamRec{..}, ref, a_sequence) -> liftIO $ do
              let good_pairs     = U.indexed             a_sequence
                  good_pairs_rev = U.indexed $ U.reverse a_sequence
                  a_fragment_type | isFirstMate  b && isPaired     b = Leading
                                  | isSecondMate b && isPaired     b = Trailing
                                  | leeHom                           = Complete
                                  | isFirstMate  b || isSecondMate b = Complete     -- old style flagging
                                  | isTrimmed    b || isMerged     b = Complete     -- new style flagging
                                  | otherwise                        = Leading

              -- basecompositon near 5' end, near 3' end
              let width  = ctx + min rng (alignedLength b_cigar `div` 2)
              mapM_ (\i -> bump (fromIntegral (ref U.!  i                   ) * maxwidth + i) acc_bc) [0 .. width-1]
              mapM_ (\i -> bump (fromIntegral (ref U.! (i + U.length ref) +6) * maxwidth + i) acc_bc) [-width .. -1]

              -- For substitutions, decide what damage class we're in:
              -- 0 - no damage, 1 - damaged 5' end, 2 - damaged 3' end, 3 - both
              let dmgbase = 2*4*4*rng * ( (if U.null a_sequence || U.head a_sequence /= (nucsC,nucsT) then 1 else 0)
                                        + (if U.null a_sequence || U.last a_sequence /= (nucsC,nucsT) then 2 else 0) )

              -- substitutions near 5' end
              let len_at_5 = case a_fragment_type of Leading  -> min rng (G.length b_seq)
                                                     Complete -> min rng (G.length b_seq `div` 2)
                                                     Trailing -> 0
              U.forM_ (U.take len_at_5 good_pairs) $
                    \(i,uv) -> withPair uv $ \j -> bump (j * rng + i + dmgbase) acc_st

              -- substitutions at CpG sites near 5' end
              U.zipWithM_
                  (\(i,(u,v)) (_,(w,z)) ->
                      when (u == nucsC && w == nucsG) $ do
                          withNs v $ \y -> bump (  y   * rng +  i ) acc_cg
                          withNs z $ \y -> bump ((y+4) * rng + i+1) acc_cg)
                  (U.take len_at_5 good_pairs) (U.drop 1 good_pairs)

              -- substitutions near 3' end
              let len_at_3 = case a_fragment_type of Leading  -> 0
                                                     Complete -> min rng (G.length b_seq `div` 2)
                                                     Trailing -> min rng (G.length b_seq)
              U.forM_ (U.take len_at_3 good_pairs_rev) $
                    \(i,uv) -> withPair uv $ \j -> bump ((17+j) * rng -i -1 + dmgbase) acc_st

              -- substitutions at CpG sites near 3' end
              U.zipWithM_
                  (\(_,(u,v)) (i,(w,z)) ->
                      when (u == nucsC && w == nucsG) $ do
                          withNs v $ \y -> bump ((y+ 9) * rng - i-2) acc_cg
                          withNs z $ \y -> bump ((y+13) * rng - i-1) acc_cg)
                  (U.drop 1 good_pairs_rev) (U.take len_at_3 good_pairs_rev)


              return ALN{..}

    let nsubsts = 4*4*rng
        mk_substs off = sequence [ (,) (n1 :-> n2) <$> U.unsafeFreeze (UM.slice ((4*i+j)*rng + off*nsubsts) rng acc_st)
                                 | (i,n1) <- zip [0..] [nucA..nucT]
                                 , (j,n2) <- zip [0..] [nucA..nucT] ]

    accs <- liftIO $ DmgStats <$> sequence [ (,) nuc <$> U.unsafeFreeze (UM.slice (i*maxwidth) maxwidth acc_bc)
                                           | (i,nuc) <- zip [2,1,3,0,4] [Just nucA, Just nucC, Just nucG, Just nucT, Nothing] ]
                              <*> sequence [ (,) nuc <$> U.unsafeFreeze (UM.slice (i*maxwidth) maxwidth acc_bc)
                                           | (i,nuc) <- zip [7,6,8,5,9] [Just nucA, Just nucC, Just nucG, Just nucT, Nothing] ]

                              <*> mk_substs 0
                              <*> mk_substs 1
                              <*> mk_substs 2
                              <*> mk_substs 3
                              <*> mk_substs 4
                              <*> mk_substs 5
                              <*> mk_substs 6
                              <*> mk_substs 7

                              <*> sequence [ (,) (n1 :-> n2) <$> U.unsafeFreeze (UM.slice ((i+j)*rng) rng acc_cg)
                                           | (i,n1) <- [(0,nucC), (4,nucG)]
                                           , (j,n2) <- zip [0..] [nucA..nucT] ]

                              <*> sequence [ (,) (n1 :-> n2) <$> U.unsafeFreeze (UM.slice ((i+j)*rng) rng acc_cg)
                                           | (i,n2) <- [(8,nucC), (12,nucG)]
                                           , (j,n1) <- zip [0..] [nucA..nucT] ]

    accs' <- accs `liftM` lift (run it')
    return $ accs' { substs5   = mconcat [ substs5 accs', substs5d5 accs', substs5d3 accs', substs5dd accs' ]
                   , substs3   = mconcat [ substs3 accs', substs3d5 accs', substs3d3 accs', substs3dd accs' ]
                   , substs5d5 = mconcat [ substs5d5 accs', substs5dd accs']
                   , substs3d5 = mconcat [ substs3d5 accs', substs3dd accs']
                   , substs5d3 = mconcat [ substs5d3 accs', substs5dd accs']
                   , substs3d3 = mconcat [ substs3d3 accs', substs3dd accs'] }
  where
    {-# INLINE withPair #-}
    withPair (Ns u, Ns v) k = case pairTab `U.unsafeIndex` fromIntegral (16*u+v) of
         j -> if j >= 0 then k j else return ()

    !pairTab = U.replicate 256 (-1) U.//
            [ (fromIntegral $ 16*u+v, x*4+y) | (Ns u,x) <- zip [nucsA, nucsC, nucsG, nucsT] [0,1,2,3]
                                             , (Ns v,y) <- zip [nucsA, nucsC, nucsG, nucsT] [0,1,2,3] ]
    {-# INLINE bump #-}
    bump i v = UM.read v i >>= UM.write v i . succ

    {-# INLINE withNs #-}
    withNs ns k | ns == nucsA = k 0
                | ns == nucsC = k 1
                | ns == nucsG = k 2
                | ns == nucsT = k 3
                | otherwise   = return ()


instance Monoid a => Monoid (DmgStats a) where
    mempty = DmgStats { basecompo5 = empty_compo
                      , basecompo3 = empty_compo
                      , substs5    = empty_subst
                      , substs3    = empty_subst
                      , substs5d5  = empty_subst
                      , substs3d5  = empty_subst
                      , substs5d3  = empty_subst
                      , substs3d3  = empty_subst
                      , substs5dd  = empty_subst
                      , substs3dd  = empty_subst
                      , substs5cpg = empty_subst
                      , substs3cpg = empty_subst
                      , stats_more = mempty }
      where
        empty_compo = [ (nuc, U.empty) | nuc <- [Just nucA, Just nucC, Just nucG, Just nucT, Nothing] ]
        empty_subst = [ (n1 :-> n2, U.empty) | n1 <- [nucA..nucT], n2 <- [nucA..nucT] ]

    a `mappend` b = DmgStats { basecompo5 = zipWith s1 (basecompo5 a) (basecompo5 b)
                             , basecompo3 = zipWith s1 (basecompo3 a) (basecompo3 b)
                             , substs5    = zipWith s2 (substs5    a) (substs5    b)
                             , substs3    = zipWith s2 (substs3    a) (substs3    b)
                             , substs5d5  = zipWith s2 (substs5d5  a) (substs5d5  b)
                             , substs3d5  = zipWith s2 (substs3d5  a) (substs3d5  b)
                             , substs5d3  = zipWith s2 (substs5d3  a) (substs5d3  b)
                             , substs3d3  = zipWith s2 (substs3d3  a) (substs3d3  b)
                             , substs5dd  = zipWith s2 (substs5dd  a) (substs5dd  b)
                             , substs3dd  = zipWith s2 (substs3dd  a) (substs3dd  b)
                             , substs5cpg = zipWith s2 (substs5cpg a) (substs5cpg b)
                             , substs3cpg = zipWith s2 (substs3cpg a) (substs3cpg b)
                             , stats_more = mappend    (stats_more a) (stats_more b) }
      where
        s1 (x, u) (z, v) | x /= z    = error "Mismatch in zip.  This is a bug."
                         | U.null u  = (x, v)
                         | U.null v  = (x, u)
                         | otherwise = (x, U.zipWith (+) u v)

        s2 (x :-> y, u) (z :-> w, v) | x /= z || y /= w = error "Mismatch in zip.  This is a bug."
                                     | U.null u         = (x :-> y, v)
                                     | U.null v         = (x :-> y, u)
                                     | otherwise        = (x :-> y, U.zipWith (+) u v)


revcom_both :: ( BamRec, U.Vector Word8, U.Vector (Nucleotides, Nucleotides) )
            -> ( BamRec, U.Vector Word8, U.Vector (Nucleotides, Nucleotides) )
revcom_both (b, ref, pairs)
    | isReversed b = ( b, revcom_ref ref, revcom_pairs pairs )
    | otherwise    = ( b,            ref,              pairs )
  where
    revcom_ref   = U.reverse . U.map (\x -> if x > 3 then x else xor x 2)
    revcom_pairs = U.reverse . U.map (compls *** compls)


-- | Reconstructs the alignment from reference, query, and cigar.  Only
-- positions where the query is not gapped are produced.
aln_from_ref :: U.Vector Word8 -> Vector_Nucs_half Nucleotides -> VS.Vector Cigar -> U.Vector NPair
aln_from_ref ref0 qry0 cig0 = U.fromList $ step ref0 qry0 cig0
  where
    step ref qry cig1
        | U.null ref || G.null qry || G.null cig1 = []
        | otherwise = case G.unsafeHead cig1 of { op :* n ->
                      case G.unsafeTail cig1 of { cig ->
                      case op of {

        Mat -> zipWith (\r q -> (nn r,q)) (G.toList (G.take n ref))
                                          (G.toList (G.take n qry)) ++ step (G.drop n ref) (G.drop n qry) cig ;
        Del ->                                                         step (G.drop n ref)           qry  cig ;
        Ins ->    map (\q -> ( gap,  q )) (G.toList (G.take n qry)) ++ step           ref  (G.drop n qry) cig ;
        SMa ->    map (\q -> ( gap,  q )) (G.toList (G.take n qry)) ++ step           ref  (G.drop n qry) cig ;
        HMa ->   replicate n (gap, nucsN)                           ++ step           ref            qry  cig ;
        Nop ->                                                         step           ref            qry  cig ;
        Pad ->                                                         step           ref            qry  cig }}}

    nn 0 = nucsT
    nn 1 = nucsC
    nn 2 = nucsA
    nn 3 = nucsG
    nn _ = nucsN


-- | Reconstructs the alignment from query, cigar, and md.  Only
-- positions where the query is not gapped are produced.
aln_from_md :: Vector_Nucs_half Nucleotides -> VS.Vector Cigar -> [MdOp] -> U.Vector NPair
aln_from_md qry0 cig0 md0 = U.fromList $ step qry0 cig0 md0
  where
    step qry cig1 md
        | G.null qry || G.null cig1 || null md = []
        | otherwise = case G.unsafeHead cig1 of op :* n -> step' qry op n (G.unsafeTail cig1) md

    step' qry  _ 0 cig             md  = step  qry      cig md
    step' qry op n cig (MdNum  0 : md) = step' qry op n cig md
    step' qry op n cig (MdDel [] : md) = step' qry op n cig md

    step' qry Mat n cig (MdNum m : md)
            | n <  m = map (\q -> (q,q)) (G.toList (G.take n qry)) ++ step  (G.drop n qry)           cig (MdNum (m-n) : md)
            | n >  m = map (\q -> (q,q)) (G.toList (G.take m qry)) ++ step' (G.drop m qry) Mat (n-m) cig                md
            | n == m = map (\q -> (q,q)) (G.toList (G.take n qry)) ++ step  (G.drop n qry)           cig                md
    step' qry Mat n cig (MdRep c : md) =         ( c, G.head qry )  : step' (G.tail   qry) Mat (n-1) cig                md
    step'   _ Mat _   _          _     = []

    step' qry Del n cig (MdDel (_:ss) : md) = step' qry Del (n-1) cig (MdDel ss : md)
    step'   _ Del _   _               _     = []

    step' qry Ins n cig                 md  = map ((,) gap) (G.toList (G.take n qry)) ++ step (G.drop n qry) cig md
    step' qry SMa n cig                 md  = map ((,) gap) (G.toList (G.take n qry)) ++ step (G.drop n qry) cig md
    step' qry HMa n cig                 md  =                replicate n (gap, nucsN) ++ step           qry  cig md
    step' qry Nop _ cig                 md  =                                            step           qry  cig md
    step' qry Pad _ cig                 md  =                                            step           qry  cig md


