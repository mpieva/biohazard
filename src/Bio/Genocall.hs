{-# LANGUAGE DeriveGeneric #-}
module Bio.Genocall where

import Bio.Adna
import Bio.Bam.Pileup
import Bio.Prelude
import Data.Aeson

import qualified Data.HashMap.Strict    as H
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

{-# INLINE simple_indel_call #-}
simple_indel_call :: (DmgToken -> Int -> Bool -> Mat44D) -> (IndelPile,IndelPile) -> (GL, [IndelVariant])
simple_indel_call get_dmg (varsF,varsR)
    | length (varsF++varsR) <= 1 = ( U.empty, [] )
    | otherwise                  = ( simple_call $ map (mkpls False) varsF ++ map (mkpls True) varsR, vars' )
  where
    vars' = IndelVariant (V_Nucs U.empty) (V_Nuc U.empty) :
            (Set.toList . Set.fromList)
                [ IndelVariant (V_Nucs $ U.fromList d)
                               (V_Nuc  $ U.fromList $ map db_call i)
                | (_q,(d,i)) <- varsF ++ varsR
                , not (null d) || not (null i) ]

    match str = zipWith $ \(DB b q dt di _) n -> let p  = get_dmg dt di str `bang` n :-> b
                                                     p' = fromQual q
                                                 in toProb $ p + p' - p * p'

    mkpls :: Bool -> (Qual, ([Nucleotides], [DamagedBase])) -> U.Vector Prob
    mkpls str (q,(d,i)) = U.fromList [ qualToProb q +
                                       if length d /= U.length dr || length i /= U.length ir
                                       then 0 else product (match str i $ U.toList ir)
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

instance ToJSON   m => ToJSON   (SubstModel_ m)
instance FromJSON m => FromJSON (SubstModel_ m)

type SubstModel = SubstModel_ Mat44D

-- | Mutable version of SubstModel, we'll probably have to accumulate in
-- this thing.
type MSubstModel = SubstModel_ MMat44D

lookupSubstModel :: SubstModel_ a -> Int -> Bool -> a
lookupSubstModel m i False
    | i >= 0 &&   i  <  V.length  (left_substs_fwd m) = V.unsafeIndex (left_substs_fwd   m)   i
    | i <  0 && (-i) <= V.length (right_substs_fwd m) = V.unsafeIndex (right_substs_fwd  m) (-i-1)
    | otherwise                                       = middle_substs_fwd m
lookupSubstModel m i True
    | i >= 0 &&   i  <  V.length  (left_substs_rev m) = V.unsafeIndex (left_substs_rev   m)   i
    | i <  0 && (-i) <= V.length (right_substs_rev m) = V.unsafeIndex (right_substs_rev  m) (-i-1)
    | otherwise                                       = middle_substs_rev m

-- Freezes a mutable substitution model into an immutable one.  Both
-- strands are combined, the result is normalized, and duplicated to
-- have a model for each strand again.
freezeSubstModel :: MSubstModel -> IO SubstModel
freezeSubstModel mm = do
    new_left   <- V.zipWithM freezeMats (left_substs_fwd   mm) (right_substs_rev  mm)
    new_middle <-            freezeMats (middle_substs_fwd mm) (middle_substs_rev mm)
    new_right  <- V.zipWithM freezeMats (right_substs_fwd  mm) (left_substs_rev   mm)

    return $ SubstModel new_left new_middle new_right
                        ( V.map complMat new_left   )
                              ( complMat new_middle )
                        ( V.map complMat new_right  )

newtype SubstModels = SubstModels (HashMap Bytes SubstModel)
  deriving (Show, Generic)

instance ToJSON SubstModels where
    toJSON (SubstModels m) = Object $ H.fromList
        [ ( decodeBytes k, toJSON v ) | (k,v) <- H.toList m ]

instance FromJSON SubstModels where
    parseJSON = withObject "map of substitution models" $ \o ->
                SubstModels . H.fromList <$> sequence
                    [ (,) (encodeBytes k) <$> parseJSON v | (k,v) <- H.toList o ]




-- | Naive SNP call; essentially the GATK model.  We compute the
-- likelihood for each base from an empirical error/damage model, then
-- hand over to 'simple_call'.  Base quality is ignored, but map quality
-- is incorporated.

{-# INLINE simple_snp_call #-}
simple_snp_call :: (DmgToken -> Int -> Bool -> Mat44D) -> (BasePile,BasePile) -> Snp_GLs
simple_snp_call get_dmg (varsF,varsR) = mk_snp_gls (simple_call $ map (mkpls False) varsF ++ map (mkpls True) varsR) ref
  where
    ref = case varsF ++ varsR of (_, DB _ _ _ _ r) : _ -> r ; _ -> nucsN
    mkpls str (qq, DB b _ dt di _) = U.generate 4 $ \n ->
                                        let x = get_dmg dt di str `bang` N (fromIntegral n) :-> b
                                        in toProb $ x + fromQual qq * (1-x)

-- | Compute @GL@ values for the simple case.  The simple case is where
-- we sample two alleles with equal probability and assume that errors
-- occur independently from each other.  This is specialized for a few
-- common cases:  two variants, because that's a typical indel; four
-- variants, because that's every SNP.

{-# INLINE simple_call #-}
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


getRow :: Int -> Mat44D -> Vec4D
getRow i (Mat44D v) = Vec4D (v U.! (4*i)) (v U.! (4*i+1)) (v U.! (4*i+2)) (v U.! (4*i+3))

setRow :: Int -> Vec4D -> Mat44D -> Mat44D
setRow i (Vec4D a b c d) (Mat44D v) = Mat44D $ v U.// [ (4*i,a), (4*i+1,b), (4*i+2,c), (4*i+3,d) ]


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
