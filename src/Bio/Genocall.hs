{-# LANGUAGE DeriveGeneric #-}
module Bio.Genocall where

import Bio.Adna
import Bio.Bam.Pileup
import Bio.Prelude

import qualified Data.Set               as Set
import qualified Data.Vector            as V
import qualified Data.Vector.Unboxed    as U

-- | Simple indel calling.  We don't bother with it too much, so here's
-- the gist:  We collect variants (simply different variants, details
-- don't matter), so \(n\) variants give rise to \((n+1)*n/2\) GL values.
-- (That's two out of \((n+1)\), the reference allele, represented here as
-- no deletion and no insertion, is there, too.)  To assign these, we
-- need a likelihood for an observed variant given an assumed genotype.
--
-- For variants of equal length, the likelihood is the sum of qualities
-- of mismatching bases, but no higher than the mapping quality.  That
-- is roughly the likelihood of getting the observed sequence even
-- though the real sequence is a different variant.  For variants of
-- different length, the likelihood is the map quality.  This
-- corresponds to the assumption that indel errors in sequencing are
-- much less likely than mapping errors.  Since this is hardly our
-- priority, the approximations are hereby declared good enough.

simple_indel_call :: IndelPile Mat44D -> (GL, [IndelVariant])
simple_indel_call   [ ] = ( U.empty, [] )
simple_indel_call   [_] = ( U.empty, [] )
simple_indel_call vars = ( simple_call $ map mkpls vars, vars' )
  where
    vars' = IndelVariant (V_Nucs U.empty) (V_Nuc U.empty) :
            (Set.toList . Set.fromList)
                [ IndelVariant (V_Nucs $ U.fromList d)
                               (V_Nuc  $ U.fromList $ map db_call i)
                | (_q,(d,i)) <- vars
                , not (null d) || not (null i) ]

    match = zipWith $ \(DB b q d _) n -> let p  = d `bang` n :-> b
                                             p' = fromQual q
                                         in toProb $ p + p' - p * p'

    mkpls :: (Qual, ([Nucleotides], [DamagedBase Mat44D])) -> U.Vector Prob
    mkpls (q,(d,i)) = U.fromList [ qualToProb q +
                                   if length d /= U.length dr || length i /= U.length ir
                                   then 0 else product (match i $ U.toList ir)
                                 | IndelVariant (V_Nucs dr) (V_Nuc ir) <- vars' ]

-- | A completely universal, completely empirical substituion model.
-- We make no attempt to distinguish damage from error.  The model is
-- cloned so we don't need to constantly flip matrices depending on
-- strand.
data SubstModel_ m = SubstModel
        { left_substs_fwd   :: {-# UNPACK #-} !(V.Vector m)
        , middle_substs_fwd ::                          !m
        , right_substs_fwd  :: {-# UNPACK #-} !(V.Vector m)
        , left_substs_rev   :: {-# UNPACK #-} !(V.Vector m)
        , middle_substs_rev ::                          !m
        , right_substs_rev  :: {-# UNPACK #-} !(V.Vector m) }
    deriving (Show, Generic)

type SubstModel = SubstModel_ Mat44D

-- | Mutable version of SubstModel, we'll probably have to accumulate in
-- this thing.
type MSubstModel = SubstModel_ MMat44D

-- | Naive SNP call; essentially the GATK model.  We compute the
-- likelihood for each base from an empirical error/damage model, then
-- hand over to 'simple_call'.  Base quality is ignored, but map quality
-- is incorporated.
--
-- XXX  Quality is no longer used, maybe it can be removed.

simple_snp_call :: BasePile Mat44D -> Snp_GLs
simple_snp_call vars = mk_snp_gls (simple_call $ map mkpls vars) ref
  where
    ref = case vars of (_, DB _ _ _ r) : _ -> r ; _ -> nucsN
    mkpls (qq, DB b _ m _) = U.generate 4 $ \n ->
                                let x = m `bang` N (fromIntegral n) :-> b
                                in toProb $ x + fromQual qq * (1-x)

-- | Compute @GL@ values for the simple case.  The simple case is where
-- we sample two alleles with equal probability and assume that errors
-- occur independently from each other.  This is specialized for a few
-- common cases:  two variants, because that's a typical indel; four
-- variants, because that's every SNP.

simple_call :: [U.Vector Prob] -> GL
simple_call [      ]                    = U.empty
simple_call (gl:gls) = case U.length gl of
    2 -> foldl' (U.zipWith (*)) (step2 gl) $ map step2 gls
              where
                step2 v = U.fromListN 3 [ x0, (x0+x1) / 2, x1 ]
                  where x0 = U.unsafeIndex v 0
                        x1 = U.unsafeIndex v 1

    4 -> foldl' (U.zipWith (*)) (step4 gl) $ map step4 gls
              where
                step4 v = U.fromListN 10 [ x0
                                         , (x0+x1)/2, x1
                                         , (x0+x2)/2, (x1+x2)/2, x2
                                         , (x0+x3)/2, (x1+x3)/2, (x2+x3)/2, x3 ]
                  where x0 = U.unsafeIndex v 0
                        x1 = U.unsafeIndex v 1
                        x2 = U.unsafeIndex v 2
                        x3 = U.unsafeIndex v 3

    _ -> foldl' (U.zipWith (*)) (step gl) $ map step gls
              where
                step !ls = U.concatMap (\i -> let hd  = U.unsafeIndex ls i
                                                  ls' = U.unsafeTake (i+1) ls
                                              in U.map (\x -> 0.5 * (hd + x)) ls'
                                       ) (U.enumFromN 0 $ U.length ls)


-- | Make a list of genotypes, each represented as a vector of allele
-- probabilities, from four possible alleles.
--
-- This makes the most sense for SNPs.  The implied order of alleles is
-- A,C,G,T, and the resulting genotype vectors can straight forwardly be
-- mutiplied with a substitution matrix to give a sensible result.
-- (Something similar for indels could be imagined, but doesn't seem all
-- that useful.  We specialize for SNPs to get simpler types and
-- efficient code.)
--
-- "For biallelic sites the ordering is: AA,AB,BB; for triallelic
-- sites the ordering is: AA,AB,BB,AC,BC,CC, etc."

mk_snp_gts :: [Vec4D]
mk_snp_gts = [ Vec4D (0.5*(a+w)) (0.5*(b+x)) (0.5*(c+y)) (0.5*(d+z))
             | as@(_:_) <- inits [ Vec4D 1 0 0 0, Vec4D 0 1 0 0, Vec4D 0 0 1 0, Vec4D 0 0 0 1 ]
             , let Vec4D a b c d = last as
             , Vec4D w x y z <- as ]

-- | SNP call according to maq/samtools/bsnp model.  The matrix k counts
-- how many errors we made, approximately.
-- XXX  Unfixable, for the time being.

{- maq_snp_call :: Int -> Double -> BasePile -> Snp_GLs
maq_snp_call theta bases = Snp_GLs (U.fromList $ map l $ mk_snp_gts ) ref
  where
    -- Bases with effective qualities in order of decreasing(!) quality.
    -- A vector based algorithm may fit here.
    bases' = sortBy (flip $ comparing db_qual)
             [ db { db_qual = mq `min` db_qual db } | (mq,db) <- bases ]

    ref = case bases of (_, DB _ _ _ r) : _ -> r ; _ -> nucsN

    nullMat = Mat44D $ U.generate 16 (const 0)
    -- L(G)
    l gt = l' gt (toProb 1) nullMat bases'

    l' !_  !acc !_ [     ] = acc
    l' !gt !acc !k (!x:xs) =
        let
            -- P(X|Q,H), a vector of four (x is fixed, h is not)
            -- this is the simple form where we set all w to 1/4
            p_x__q_h_ :: Vec4D
            p_x__q_h_ = vecNucs $ \h -> 0.25 * fromQualRaised (theta ** (k `bang` h :-> db_call x)) (db_qual x)

            -- eh, this is cumbersome... what was I thinking?!
            p_x__q_h  :: Vec4D
            p_x__q_h  = vecZipNucs (\p h -> if db_call x == h then 1 + p - vecSum p_x__q_h_ else p) p_x__q_h_

            -- P(H|X), again a vector of four
            p_x__q    = dot p_x__q_h dg
            p_h__x    = vecZip (\p p_h -> p / p_x__q * p_h) p_x__q_h dg
            dg        = db_dmg x `multmv` gt

            kk = vecZip (+) (getRow (fromIntegral . unN $ db_call x) k) p_h__x
            k' = setRow (fromIntegral . unN $ db_call x) kk k

            acc' = acc * toProb p_x__q
            meh = vecNucs $ \h -> k `bang` h :-> db_call x
        in {- trace (unlines ["gt " ++ show gt
                          ,"p(x|q,h) " ++ show p_x__q_h
                          ,"dg " ++ show dg ++ ", call = " ++ show (db_call x)
                          ,"p(h|x) " ++ show p_h__x
                          ,"k  " ++ show k
                          ,"k' " ++ show k'
                          ,"meh " ++ show meh]) $ -} l' gt acc' k' xs -}

getRow :: Int -> Mat44D -> Vec4D
getRow i (Mat44D v) = Vec4D (v U.! (4*i)) (v U.! (4*i+1)) (v U.! (4*i+2)) (v U.! (4*i+3))

setRow :: Int -> Vec4D -> Mat44D -> Mat44D
setRow i (Vec4D a b c d) (Mat44D v) = Mat44D $ v U.// [ (4*i,a), (4*i+1,b), (4*i+2,c), (4*i+3,d) ]


--  Error model with dependency parameter.  Since both strands are
-- supposed to still be independent, we feed in only one pile, and
-- later combine both calls.  XXX What's that doing HERE?!

type Calls = Pile' Snp_GLs (GL, [IndelVariant])

-- | This pairs up GL values and the reference allele.  When
-- constructing it, we make sure the GL values are in the correct order
-- if the reference allele is listed first.
data Snp_GLs = Snp_GLs { snp_gls :: !GL, snp_refbase :: !Nucleotides }
    deriving Show

mk_snp_gls :: GL -> Nucleotides -> Snp_GLs
mk_snp_gls gl ref | U.length gl /= 10 = error "only diploid genomes are supported!"
                  | otherwise         = Snp_GLs gl ref

data Vec4D = Vec4D {-# UNPACK #-} !Double {-# UNPACK #-} !Double {-# UNPACK #-} !Double {-# UNPACK #-} !Double

vecNucs :: (Nucleotide -> Double) -> Vec4D
vecNucs f = Vec4D (f nucA) (f nucC) (f nucG) (f nucT)

vecSum :: Vec4D -> Double
vecSum (Vec4D a b c d)  = a + b + c + d

dot :: Vec4D -> Vec4D -> Double
dot (Vec4D a b c d) (Vec4D w x y z) = a*w + b*x + c*y + d*z

multmv :: Mat44D -> Vec4D -> Vec4D
multmv m v = Vec4D (dot (getRow 0 m) v) (dot (getRow 1 m) v)
                   (dot (getRow 2 m) v) (dot (getRow 3 m) v)

vecZip :: (Double -> Double -> Double) -> Vec4D -> Vec4D -> Vec4D
vecZip f (Vec4D a b c d) (Vec4D w x y z) = Vec4D (f a w) (f b x) (f c y) (f d z)

vecZipNucs :: (Double -> Nucleotide -> Double) -> Vec4D -> Vec4D
vecZipNucs f (Vec4D a b c d) = Vec4D (f a nucA) (f b nucC) (f c nucG) (f d nucT)
