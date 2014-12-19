{-# LANGUAGE BangPatterns, RecordWildCards #-}
module Bio.Genocall.Adna where

import Bio.Base
import Bio.Bam ( BamRaw, br_isReversed, br_l_seq )
import Data.Vec hiding ( map )

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


-- | A 'DamageModel' is a function that gives substitution matrices for
-- each position in a read.  Its application yields a sequence of
-- substitution matrices exactly as long a the read itself.  Though
-- typically not done, the model can take the sequence into account.
-- It will usually be dependendent on the position within the read,
-- though.

type DamageModel a = BamRaw -> [Mat44 a]

data To = Nucleotide :-> Nucleotide

infix 9 :->
infix 8 !

-- | Convenience function to access a substitution matrix that has a
-- mnemonic reading.
(!) :: Mat44D -> To -> Double
(!) m (N x :-> N y) = getElem (fromIntegral x) $ getElem (fromIntegral y) m

-- | 'DamageModel' for undamaged DNA.  The likelihoods follow directly
-- from the quality score.  This needs elaboration to see what to do
-- with amibiguity codes (even though those haven't actually been
-- observed in the wild).

{-# SPECIALIZE noDamage :: DamageModel Double #-}
noDamage :: Num a => DamageModel a
noDamage r = replicate (br_l_seq r) identity


-- | 'DamageModel' for single stranded library prep.  Only one kind of
-- damage occurs (C to T), it occurs at low frequency ('delta_ds')
-- everywhere, at high frequency ('delta_ss') in single stranded parts,
-- and the overhang length is distributed exponentially with parameter
-- 'lambda' at the 5' end and 'kappa' at the 3' end.

data SsDamageParameters float = SSD { ssd_sigma  :: !float         -- deamination rate in ss DNA
                                    , ssd_delta  :: !float         -- deamination rate in ds DNA
                                    , ssd_lambda :: !float         -- expected overhang length at 5' end
                                    , ssd_kappa  :: !float }       -- expected overhang length at 3' end
  deriving Show

-- Forward strand first, C->T only; reverse strand next, G->A instead

{-# SPECIALIZE ssDamage :: SsDamageParameters Double -> DamageModel Double #-}
ssDamage :: Fractional a => SsDamageParameters a -> DamageModel a
ssDamage SSD{..} r = let mat = if br_isReversed r then ssd_rev else ssd_fwd
                     in map mat [0 .. br_l_seq r-1]
  where
    prob5 = abs ssd_lambda / (1 + abs ssd_lambda)
    prob3 = abs ssd_kappa  / (1 + abs ssd_kappa)

    ssd_fwd i = vec4 ( vec4  1   0   0   0 )
                     ( vec4  0 (1-p) 0   0 )
                     ( vec4  0   0   1   0 )
                     ( vec4  0   p   0   1 )
      where
        !lam5 = prob5 ^ (1+i)
        !lam3 = prob3 ^ (br_l_seq r - i)
        !lam  = lam3 + lam5 - lam3 * lam5
        !p    = ssd_sigma * lam + ssd_delta * (1-lam)

    ssd_rev i = vec4 ( vec4  1   0   p   0 )
                     ( vec4  0   1   0   0 )
                     ( vec4  0   0 (1-p) 0 )
                     ( vec4  0   0   0   1 )
      where
        !lam5 = prob5 ^ (br_l_seq r - i)
        !lam3 = prob3 ^ (1+i)
        !lam  = lam3 + lam5 - lam3 * lam5
        !p    = ssd_sigma * lam + ssd_delta * (1-lam)


data DsDamageParameters float = DSD { dsd_sigma  :: !float         -- deamination rate in ss DNA
                                    , dsd_delta  :: !float         -- deamination rate in ds DNA
                                    , dsd_lambda :: !float }       -- expected overhang length
  deriving Show

-- | 'DamageModel' for double stranded library.  We get C->T damage at
-- the 5' end and G->A at the 3' end.  Everything is symmetric, and
-- therefore the orientation of the aligned read doesn't matter either.
--
-- Parameterization is stolen from @mapDamage 2.0@, for the most part.
-- We have a deamination rate each for ss and ds dna and average
-- overhang length parameters for both ends.  (Without UDG treatment,
-- those will be equal.  With UDG, those are much smaller and in fact
-- don't literally represent overhangs.)
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

{-# SPECIALIZE dsDamage :: DsDamageParameters Double -> DamageModel Double #-}
dsDamage :: Fractional a => DsDamageParameters a -> DamageModel a
dsDamage DSD{..} r = map mat [0 .. br_l_seq r-1]
  where
    prob = abs dsd_lambda / (1 + abs dsd_lambda)

    mat i = vec4 ( vec4  1     0     q     0 )
                 ( vec4  0   (1-p)   0     0 )
                 ( vec4  0     0   (1-q)   0 )
                 ( vec4  0     p     0     1 )
      where
        p    = dsd_sigma * lam5 + dsd_delta * (1-lam5)
        q    = dsd_sigma * lam3 + dsd_delta * (1-lam3)
        lam5 = prob ^ (         1 + i)
        lam3 = prob ^ (br_l_seq r - i)

{-# INLINE vec4 #-}
vec4 :: a -> a -> a -> a -> Vec4 a
vec4 a b c d = a :. b :. c :. d :. ()
