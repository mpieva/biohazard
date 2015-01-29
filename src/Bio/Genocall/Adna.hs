{-# LANGUAGE BangPatterns, RecordWildCards #-}
module Bio.Genocall.Adna where

import Bio.Base
import Data.Vec
import qualified Data.Vector as V

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
-- each position in a read.  The 'DamageModel' can depend on the length
-- of the read and whether its alignment is reversed.  In practice, we
-- should probably memoize precomputed damage models somehow.

type DamageModel a = Bool -> Int -> V.Vector (Mat44 a)

data To = Nucleotide :-> Nucleotide

infix 9 :->
infix 8 !

-- | Convenience function to access a substitution matrix that has a
-- mnemonic reading.
{-# INLINE (!) #-}
(!) :: Mat44D -> To -> Double
(!) m (N x :-> N y) = getElem (fromIntegral x) $ getElem (fromIntegral y) m

-- | 'DamageModel' for undamaged DNA.  The likelihoods follow directly
-- from the quality score.  This needs elaboration to see what to do
-- with amibiguity codes (even though those haven't actually been
-- observed in the wild).

{-# SPECIALIZE noDamage :: DamageModel Double #-}
noDamage :: Num a => DamageModel a
noDamage _ l = V.replicate l identity


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
ssDamage SSD{..} r l = V.generate l $ if r then ssd_rev else ssd_fwd
  where
    prob5 = abs ssd_lambda / (1 + abs ssd_lambda)
    prob3 = abs ssd_kappa  / (1 + abs ssd_kappa)

    ssd_fwd i = vec4 ( vec4  1   0   0   0 )
                     ( vec4  0 (1-p) 0   0 )
                     ( vec4  0   0   1   0 )
                     ( vec4  0   p   0   1 )
      where
        !lam5 = prob5 ^ (1+i)
        !lam3 = prob3 ^ (l-i)
        !lam  = lam3 + lam5 - lam3 * lam5
        !p    = ssd_sigma * lam + ssd_delta * (1-lam)

    ssd_rev i = vec4 ( vec4  1   0   p   0 )
                     ( vec4  0   1   0   0 )
                     ( vec4  0   0 (1-p) 0 )
                     ( vec4  0   0   0   1 )
      where
        !lam5 = prob5 ^ (l-i)
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
dsDamage DSD{..} _ l = V.generate l mat
  where
    prob = abs dsd_lambda / (1 + abs dsd_lambda)

    mat i = vec4 ( vec4  1     0     q     0 )
                 ( vec4  0   (1-p)   0     0 )
                 ( vec4  0     0   (1-q)   0 )
                 ( vec4  0     p     0     1 )
      where
        p    = dsd_sigma * lam5 + dsd_delta * (1-lam5)
        q    = dsd_sigma * lam3 + dsd_delta * (1-lam3)
        lam5 = prob ^ (1+i)
        lam3 = prob ^ (l-i)

{-# INLINE vec4 #-}
vec4 :: a -> a -> a -> a -> Vec4 a
vec4 a b c d = a :. b :. c :. d :. ()

memoDamageModel :: DamageModel a -> DamageModel a
memoDamageModel f = \r l -> if l > 512 || l < 0 then f r l
                            else if r then V.unsafeIndex rev l
                            else           V.unsafeIndex fwd l
  where
    rev = V.generate 512 $ f True
    fwd = V.generate 512 $ f False

