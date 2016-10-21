module Bio.Genocall where

import Bio.Adna
import Bio.Bam.Pileup
import Bio.Prelude

import qualified Data.Set               as Set
import qualified Data.Vector.Unboxed    as V

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

simple_indel_call :: Int -> IndelPile -> (GL, [IndelVariant])
simple_indel_call      _  [ ] = ( V.empty, [] )
simple_indel_call      _  [_] = ( V.empty, [] )
simple_indel_call ploidy vars = ( simple_call ploidy $ map mkpls vars, vars' )
  where
    vars' = IndelVariant (V_Nucs V.empty) (V_Nuc V.empty) :
            (Set.toList . Set.fromList)
                [ IndelVariant (V_Nucs $ V.fromList d)
                               (V_Nuc  $ V.fromList $ map db_call i)
                | (_q,(d,i)) <- vars
                , not (null d) || not (null i) ]

    match = zipWith $ \(DB b q _ m) n -> let p  = m `bang` n :-> b
                                             p' = fromQual q
                                         in toProb $ p + p' - p * p'

    mkpls :: (Qual, ([Nucleotides], [DamagedBase])) -> V.Vector Prob
    mkpls (q,(d,i)) = V.fromList [ qualToProb q +
                                   if length d /= V.length dr || length i /= V.length ir
                                   then 0 else product (match i $ V.toList ir)
                                 | IndelVariant (V_Nucs dr) (V_Nuc ir) <- vars' ]

-- | Naive SNP call; essentially the GATK model.  We create a function
-- that computes a likelihood for a given base, then hand over to simple
-- call.  Since everything is so straight forward, this works even in
-- the face of damage.

simple_snp_call :: (Qual -> Double) -> Int -> BasePile -> Snp_GLs
simple_snp_call from_qual ploidy vars = snp_gls (simple_call ploidy $ map mkpls vars) ref
  where
    ref = case vars of (_, DB _ _ r _) : _ -> r ; _ -> nucsN
    mkpls (q, DB b qq _ m) = V.generate 4 $ \n ->
                                let x = m `bang` N (fromIntegral n) :-> b
                                in toProb $ x + pe*(s-x)
      where
        !p1 = from_qual q
        !p2 = from_qual qq
        !pe = p1 + p2 - p1*p2
        !s  = sum [ m `bang` N n :-> b | n <- [0..3] ] / 4

-- | Compute @GL@ values for the simple case.  The simple case is where
-- we sample 'ploidy' alleles with equal probability and assume that
-- errors occur independently from each other.
--
-- XXX  This eats up ~40% of total runtime; it *screams out* for
-- specialization to diploidy and four alleles (common SNPs) and maybe
-- diploidy and two alleles (common indels).

simple_call :: Int -> [V.Vector Prob] -> GL
simple_call     !_ [      ]                    = V.empty
simple_call      0        _                    = V.empty
simple_call      1 (gl:gls)                    = foldl' (V.zipWith (*)) gl gls
simple_call      2 (gl:gls) | V.length gl == 2 = foldl' (V.zipWith (*)) (step2 gl) $ map step2 gls
  where
    step2 v = V.fromListN 3 [ x0, (x0+x1) / 2, x1 ]
      where x0 = V.unsafeIndex v 0
            x1 = V.unsafeIndex v 1

simple_call      2 (gl:gls) | V.length gl == 4 = foldl' (V.zipWith (*)) (step4 gl) $ map step4 gls
  where
    step4 v = V.fromListN 10 [ x0
                             , (x0+x1)/2, x1
                             , (x0+x2)/2, (x1+x2)/2, x2
                             , (x0+x3)/2, (x1+x3)/2, (x2+x3)/2, x3 ]
      where x0 = V.unsafeIndex v 0
            x1 = V.unsafeIndex v 1
            x2 = V.unsafeIndex v 2
            x3 = V.unsafeIndex v 3

simple_call ploidy (gl:gls)                    = foldl' (V.zipWith (*)) (step gl) $ map step gls
  where
    !mag = recip $ toProb (fromIntegral ploidy)

    step = mk_pls ploidy 0 . V.map (mag *)

    -- Meh.  Pointless, but happens to be the unit.
    mk_pls 0 !acc  !_ = V.singleton acc

    -- Okay, we sample ONE allele.  Likelihood of the data is simply the
    -- GL value that was passed to us.
    mk_pls 1 !acc !ls = V.map (acc +) ls

    -- We extend the genotype and sample another allele.

    mk_pls n !acc !ls = V.concatMap (\i ->
                                 let hd  = V.unsafeIndex ls i
                                     ls' = V.unsafeTake (i+1) ls
                                 in (mk_pls (n-1) (hd + acc) ls'))
                             (V.enumFromN 0 l)
            where l = V.length ls


-- | Make a list of genotypes, each represented as a vector of allele
-- probabilities, from ploidy and four possible alleles.
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

mk_snp_gts :: Int -> [Vec4D]
mk_snp_gts ploidy = go ploidy alleles
  where
    !mag = recip $ fromIntegral ploidy
    alleles = [ Vec4D 1 0 0 0, Vec4D 0 1 0 0, Vec4D 0 0 1 0, Vec4D 0 0 0 1 ]

    -- 'go p' as returns all p-ploid genotypes that can be made from the
    -- alleles 'as', in the order in which they appear in VCF.
    -- So, that's
    --   - all (p-1)-ploid genotypes that can be made from 1 allele, plus allele 0        (AA)
    --   - all (p-1)-ploid genotypes that can be made from 2 alleles, plus allele 1       (AC,CC)
    --     ...
    --
    --   - there's one 0-ploid genotype: the zero vector
    --   - the genotypes that can be made from 0 alleles is an empty list

    go !p as | p == 0    = [ Vec4D 0 0 0 0 ]
             | otherwise = [ macsvv mag (last as') gt {-gt + mag * last as'-} | as'@(_:_) <- inits as, gt <- go (p-1) as' ]

-- | SNP call according to maq/samtools/bsnp model.  The matrix k counts
-- how many errors we made, approximately.

maq_snp_call :: Int -> Double -> BasePile -> Snp_GLs
maq_snp_call ploidy theta bases = snp_gls (V.fromList $ map l $ mk_snp_gts ploidy) ref
  where
    -- Bases with effective qualities in order of decreasing(!) quality.
    -- A vector based algorithm may fit here.
    bases' = sortBy (flip $ comparing db_qual)
             [ db { db_qual = mq `min` db_qual db } | (mq,db) <- bases ]

    ref = case bases of (_, DB _ _ r _) : _ -> r ; _ -> nucsN

    nullMat = Mat44D $ V.generate 16 (const 0)
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
                          ,"meh " ++ show meh]) $ -} l' gt acc' k' xs

getRow :: Int -> Mat44D -> Vec4D
getRow i (Mat44D v) = Vec4D (v V.! (4*i)) (v V.! (4*i+1)) (v V.! (4*i+2)) (v V.! (4*i+3))

setRow :: Int -> Vec4D -> Mat44D -> Mat44D
setRow i (Vec4D a b c d) (Mat44D v) = Mat44D $ v V.// [ (4*i,a), (4*i+1,b), (4*i+2,c), (4*i+3,d) ]


{-
smoke_test :: IO ()
smoke_test =
    -- decodeAnyBamFile "/mnt/datengrab/test.bam" >=> run $ \_hdr ->
    -- enumPure1Chunk crap_data >=> run $
    -- joinI $ filterStream ((/=) (Q 0) . br_mapq) $
    -- joinI $ pileup (dsDamage $ DSD 0.9 0.02 0.3) $ -- noDamage $
    joinI $ pileup (ssDamage $ SSD 0.9 0.02 0.3 0.5) $ -- noDamage $
    -- joinI $ takeStream 5 $ mapStreamM_ print
    -- joinI $ filterStream ((> 0) . either vc_mapq0 vc_mapq0) $
    joinI $ takeStream 5000 $ mapStreamM_ call_and_print
  where
    call_and_print (Right ic) = put . showCall show_indels . fmap (simple_indel_call 2) $ ic
    call_and_print (Left  bc) = put . showCall show_bases  . fmap (simple_snp_call   2) $ bc

    put f = putStr $ f "\n"

    show_bases :: () -> ShowS
    show_bases () = (++) "A,C,G,T"

    show_indels :: IndelVars -> ShowS
    show_indels = (++) . intercalate "," . map show_indel

    show_indel :: (Int, [Nucleotide]) -> String
    show_indel (d, ins) = shows ins $ '-' : show d
-}

--  Error model with dependency parameter.  Since both strands are
-- supposed to still be independent, we feed in only one pile, and
-- later combine both calls.  XXX What's that doing HERE?!

type Calls = Pile' Snp_GLs (GL, [IndelVariant])

-- | This pairs up GL values and the reference allele.  When
-- constructing it, we make sure the GL values are in the correct order
-- if the reference allele is listed first.
data Snp_GLs = Snp_GLs !GL !Nucleotides
    deriving Show

snp_gls :: GL -> Nucleotides -> Snp_GLs
snp_gls pls ref | ref == nucsT = Snp_GLs (pls `V.backpermute` V.fromList [9,6,0,7,1,2,8,3,4,5]) ref
                | ref == nucsG = Snp_GLs (pls `V.backpermute` V.fromList [5,3,0,4,1,2,8,6,7,9]) ref
                | ref == nucsC = Snp_GLs (pls `V.backpermute` V.fromList [2,1,0,4,3,5,7,6,8,9]) ref
                | otherwise    = Snp_GLs pls ref


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

-- | 'macsvv s u v' computes s*u+v.
macsvv :: Double -> Vec4D -> Vec4D -> Vec4D
macsvv s (Vec4D a b c d) (Vec4D w x y z) = Vec4D (s*a+w) (s*b+x) (s*c+y) (s*d+z)

vecZip :: (Double -> Double -> Double) -> Vec4D -> Vec4D -> Vec4D
vecZip f (Vec4D a b c d) (Vec4D w x y z) = Vec4D (f a w) (f b x) (f c y) (f d z)

vecZipNucs :: (Double -> Nucleotide -> Double) -> Vec4D -> Vec4D
vecZipNucs f (Vec4D a b c d) = Vec4D (f a nucA) (f b nucC) (f c nucG) (f d nucT)
