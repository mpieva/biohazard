{-# LANGUAGE BangPatterns, RecordWildCards #-}
module Bio.Adna where

import Bio.Base
import Bio.Bam ( BamRaw, br_isReversed, br_l_seq )
import Numeric ( showFFloat )

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

-- | For likelihoods for bases @A, C, G, T@, one quality score.  Note
-- that we cannot roll the quality score into the probabilities:  the
-- @DB@ only describes changes that happened before sequencing, then
-- quality describes those that happen while sequencing.  The latter
-- behave differently (notably, they repeat semi-systematically).
data DamagedBase = DB !Double !Double !Double !Double !Qual

instance Show DamagedBase where
    showsPrec _ (DB a c g t q) = unwordS $ [ showFFloat (Just 2) p | p <- [a,c,g,t] ] ++ [ shows q ]
      where unwordS = foldr1 (\u v -> u . (:) ' ' . v)


-- | A 'DamageModel' is a function that gives likelihoods for all
-- possible four bases, given the sequenced base, the quality, and the
-- position in the read.  That means it can't take the actual alignment
-- into account... but nobody seems too keen on doing that anyway.
type DamageModel = BamRaw           -- ^ the read itself
                -> Int              -- ^ position in read
                -> Nucleotide       -- ^ base
                -> Qual             -- ^ quality score
                -> DamagedBase      -- ^ results in four likelihoods

-- | 'DamageModel' for undamaged DNA.  The likelihoods follow directly
-- from the quality score.  This needs elaboration to see what to do
-- with amibiguity codes (even though those haven't actually been
-- observed in the wild).
noDamage :: DamageModel
noDamage _ _ b | b == nucA = DB 1 0 0 0
               | b == nucC = DB 0 1 0 0
               | b == nucG = DB 0 0 1 0
               | b == nucT = DB 0 0 0 1
               | otherwise = DB f f f f where f = 0.25


-- | 'DamageModel' for single stranded library prep.  Only one kind of
-- damage occurs (C to T), it occurs at low frequency ('delta_ds')
-- everywhere, at high frequency ('delta_ss') in single stranded parts,
-- and the overhang length is distributed exponentially with parameter
-- 'lambda'.

data SsDamageParameters = SSD { ssd_sigma  :: !Double         -- deamination rate in ss DNA
                              , ssd_delta  :: !Double         -- deamination rate in ds DNA
                              , ssd_lambda :: !Double         -- expected overhang length at 5' end
                              , ssd_kappa  :: !Double }       -- expected overhang length at 3' end
  deriving Show

-- Forward strand first, C->T only; reverse strand next, G->A instead
-- N averages over all others, horizontally(!).  Distance from end
-- depends on strand, too :(
ssDamage :: SsDamageParameters -> DamageModel
ssDamage SSD{..} r i b = if br_isReversed r then ssd_rev else ssd_fwd
  where
    dq x y z w = DB (0.25*x) (0.25*y) (0.25*z) (0.25*w)
    prob5 = abs ssd_lambda / (1 + abs ssd_lambda)
    prob3 = abs ssd_kappa  / (1 + abs ssd_kappa)

    ssd_fwd | b == nucA = DB   1   0   0   0
            | b == nucC = DB   0 (1-p) 0   0
            | b == nucG = DB   0   0   1   0
            | b == nucT = DB   0   p   0   1
            | otherwise = dq   1 (1-p) 1 (1+p)
      where
        !lam5 = prob5 ^ (1+i)
        !lam3 = prob3 ^ (br_l_seq r - i)
        !lam  = lam3 + lam5 - lam3 * lam5
        !p    = ssd_sigma * lam + ssd_delta * (1-lam)

    ssd_rev | b == nucA = DB   1   0   p   0
            | b == nucC = DB   0   1   0   0
            | b == nucG = DB   0   0 (1-p) 0
            | b == nucT = DB   0   0   0   1
            | otherwise = dq (1+p) 1 (1-p) 1
      where
        !lam5 = prob5 ^ (br_l_seq r - i)
        !lam3 = prob3 ^ (1+i)
        !lam  = lam3 + lam5 - lam3 * lam5
        !p    = ssd_sigma * lam + ssd_delta * (1-lam)


data DsDamageParameters = DSD { dsd_sigma  :: !Double         -- deamination rate in ss DNA
                              , dsd_delta  :: !Double         -- deamination rate in ds DNA
                              , dsd_lambda :: !Double }       -- expected overhang length
  deriving Show

-- | 'DamageModel' for double stranded library.  We get C->T damage at
-- the 5' end and G->A at the 3' end.  Everything is symmetric, and
-- therefore the orientation of the aligned read doesn't matter either.
--
-- Parameterization is stolen from @mapDamage 2.0@, for the most part.
-- We have a deamination rate each for ss and ds dna and average
-- overhang length parameters for both ends.  (Without UDG treatment,
-- those will be equal.  With UDG, those are much smaller and unequal,
-- and in fact don't literally represent overhangs.)
--
-- L({A,C,G,T}|A) = {1,0,0,0}
-- L({A,C,G,T}|C) = {0,1-p,0,p}
-- L({A,C,G,T}|G) = {0,0,1,0}
-- L({A,C,G,T}|T) = {0,0,0,1}
--
-- Here, p is the probability of deamination, which is dsd_delta if
-- double strande, dsd_sigma if single stranded.  The probability of
-- begin single stranded is prob ^ (i+1).  This gives an average
-- overhang length of (prob / (1-prob)).  We invert this and define
-- dsd_lambda as the expected overhang length.

dsDamage :: DsDamageParameters -> DamageModel
dsDamage DSD{..} r i b | b == nucA = DB    1     0     q     0
                       | b == nucC = DB    0   (1-p)   0     0
                       | b == nucG = DB    0     0   (1-q)   0
                       | b == nucT = DB    0     p     0     1
                       | otherwise = dbq (1+q) (1-p) (1-q) (1+p)
  where
    prob = abs dsd_lambda / (1 + abs dsd_lambda)

    p    = dsd_sigma * lam5 + dsd_delta * (1-lam5)
    q    = dsd_sigma * lam3 + dsd_delta * (1-lam3)
    lam5 = prob ^ (         1 + i)
    lam3 = prob ^ (br_l_seq r - i)
    dbq x y z w = DB (0.25*x) (0.25*y) (0.25*z) (0.25*w)

