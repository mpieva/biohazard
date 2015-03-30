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
  deriving (Read, Show)

-- | Generic substitution matrix, has C->T and G->A deamination as
-- parameters.  Setting 'p' or 'q' to 0 as appropriate makes this apply
-- to the single stranded or undamaged case.

genSubstMat :: Fractional a => a -> a -> Mat44 a
genSubstMat p q = vec4 ( vec4  1   0     q   0 )
                       ( vec4  0 (1-p)   0   0 )
                       ( vec4  0   0   (1-q) 0 )
                       ( vec4  0   p     0   1 )

-- Forward strand first, C->T only; reverse strand next, G->A instead

{-
{-# SPECIALIZE ssDamage :: SsDamageParameters Double -> DamageModel Double #-}
ssDamage :: Fractional a => SsDamageParameters a -> DamageModel a
ssDamage SSD{..} r l = V.generate l $ if r then ssd_rev else ssd_fwd
  where
    ssd_fwd i = genSubstMat p 0
      where
        !lam5 = ssd_lambda ^ (1+i)
        !lam3 = ssd_kappa ^ (l-i)
        !lam  = lam3 + lam5 - lam3 * lam5
        !p    = ssd_sigma * lam + ssd_delta * (1-lam)

    ssd_rev i = genSubstMat 0 p
      where
        !lam5 = ssd_lambda ^ (l-i)
        !lam3 = ssd_kappa ^ (1+i)
        !lam  = lam3 + lam5 - lam3 * lam5
        !p    = ssd_sigma * lam + ssd_delta * (1-lam)



{-# SPECIALIZE dsDamage :: DsDamageParameters Double -> DamageModel Double #-}
dsDamage :: Fractional a => DsDamageParameters a -> DamageModel a
dsDamage DSD{..} _ l = V.generate l mat
  where
    mat i = genSubstMat p q
      where
        p    = dsd_sigma * lam5 + dsd_delta * (1-lam5)
        q    = dsd_sigma * lam3 + dsd_delta * (1-lam3)
        lam5 = dsd_lambda ^ (1+i)
        lam3 = dsd_lambda ^ (l-i)
-}

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

{-# SPECIALIZE univDamage :: DamageParameters Double -> DamageModel Double #-}
univDamage :: Fractional a => DamageParameters a -> DamageModel a
univDamage DP{..} r l = V.generate l mat
  where
    mat i = genSubstMat (p1+p2) (q1+q2)
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

