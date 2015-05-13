{-# LANGUAGE BangPatterns #-}
module Bio.Genocall where

import Debug.Trace

import Bio.Bam.Pileup
import Bio.Base
import Bio.Genocall.Adna
import Control.Applicative
import Data.Foldable hiding ( sum, product )
import Data.List ( inits, tails, sortBy )
import Data.Ord
import Data.Vec.Base ( (:.)(..) )
import Data.Vec.LinAlg
import Data.Vec.Packed

import qualified Data.Set               as Set
import qualified Data.Vector.Unboxed    as V
import qualified Data.Vec               as Vec

-- | Simple indel calling.  We don't bother with it too much, so here's
-- the gist:  We collect variants (simply different variants, details
-- don't matter), so @n@ variants give rise to (n+1)*n/2 GL values.
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

simple_indel_call :: Int -> IndelPile -> (GL, [IndelVariant])
simple_indel_call _ [ ] = (V.empty, [])
simple_indel_call _ [_] = (V.empty, [])
simple_indel_call ploidy vars = (simple_call ploidy mkpls vars, vars')
  where
    vars' = Set.toList $ Set.fromList
            [ IndelVariant d (V_Nuc $ V.fromList $ map db_call i) | (_q,(d,i)) <- vars ]

    match = zipWith $ \(DB b q m) n -> let p  = m ! n :-> b
                                           p' = fromQual q
                                       in toProb $ p + p' - p * p'

    mkpls (q,(d,i)) = let !q' = qualToProb q
                      in [ if d /= dr || length i /= V.length ir
                           then q' else q' + product (match i $ V.toList ir)
                         | IndelVariant dr (V_Nuc ir) <- vars' ]

-- | Naive SNP call; essentially the GATK model.  We create a function
-- that computes a likelihood for a given base, then hand over to simple
-- call.  Since everything is so straight forward, this works even in
-- the face of damage.

simple_snp_call :: Int -> BasePile -> GL
simple_snp_call ploidy vars = simple_call ploidy mkpls vars
  where
    mkpls (q, DB b qq m) = [ toProb $ x + pe*(s-x) | n <- [0..3], let x = m ! N n :-> b ]
      where
        !p1 = fromQual q
        !p2 = fromQual qq
        !pe = p1 + p2 - p1*p2
        !s  = sum [ m ! N n :-> b | n <- [0..3] ] / 4

-- | Compute @GL@ values for the simple case.  The simple case is where
-- we sample 'ploidy' alleles with equal probability and assume that
-- errors occur independently from each other.
--
-- The argument 'pls' is a function that computes the likelihood for
-- getting the current read, for every variant assuming that variant was
-- sampled.
--
-- NOTE, this may warrant specialization to diploidy and four alleles
-- (common SNPs) and diploidy and two alleles (common indels).

simple_call :: Int -> (a -> [Prob Double]) -> [a] -> GL
simple_call ploidy pls = foldl1' (V.zipWith (*)) . map step
  where
    foldl1' _ [    ] = V.singleton 1
    foldl1' f (a:as) = foldl' f a as

    !mag = recip $ toProb (fromIntegral ploidy)

    -- XXX This could probably be simplified given the mk_pls function
    -- below.
    step = V.fromList . map (* mag) . reverse . mk_pls ploidy . reverse . pls

    -- Meh.  Pointless, but happens to be the unit.
    mk_pls 0  _ = return 0

    -- Okay, we sample ONE allele.  Likelihood of the data is simply the
    -- GL value that was passed to us.
    mk_pls 1 ls = ls

    -- We extend the genotype and sample another allele.
    mk_pls n ls = do ls'@(hd:_) <- tails ls
                     (+) hd <$> mk_pls (n-1) ls'


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
             | otherwise = [ gt + mag * last as' | as'@(_:_) <- inits as, gt <- go (p-1) as' ]

-- | SNP call according to maq/samtools/bsnp model.  The matrix k counts
-- how many errors we made, approximately.

maq_snp_call :: Int -> Double -> BasePile -> GL
maq_snp_call ploidy theta bases = V.fromList $ map l $ mk_snp_gts ploidy
  where
    -- Bases with effective qualities in order of decreasing(!) quality.
    -- A vector based algorithm may fit here.
    bases' = sortBy (flip $ comparing db_qual)
             [ db { db_qual = mq `min` db_qual db } | (mq,db) <- bases ]

    everynuc :: Vec.Vec4 Nucleotide
    everynuc = nucA :. nucC :. nucG :. nucT :. ()

    -- L(G)
    l gt = l' gt (toProb 1) (0 :: Mat44D) bases'

    l' !_  !acc !_ [     ] = acc
    l' !gt !acc !k (!x:xs) =
        let
            -- P(X|Q,H), a vector of four (x is fixed, h is not)
            -- this is the simple form where we set all w to 1/4
            p_x__q_h_ = Vec.map (\h -> 0.25 * fromQualRaised (theta ** (k ! h :-> db_call x)) (db_qual x)) everynuc

            -- eh, this is cumbersome... what was I thinking?!
            p_x__q_h  = Vec.zipWith (\p h -> if db_call x == h then 1 + p - Vec.sum p_x__q_h_ else p) p_x__q_h_ everynuc

            -- P(H|X), again a vector of four
            p_x__q   = dot p_x__q_h dg
            p_h__x   = Vec.zipWith (\p p_h -> p / p_x__q * p_h) p_x__q_h dg
            dg = db_dmg x `multmv` gt

            kk = Vec.getElem (fromIntegral . unN $ db_call x) k + pack p_h__x
            k' = Vec.setElem (fromIntegral . unN $ db_call x) kk k

            acc' = acc * toProb p_x__q
            meh = Vec.map (\h -> k ! h :-> db_call x) everynuc -- XXX
        in {- trace (unlines ["gt " ++ show gt
                          ,"p(x|q,h) " ++ show p_x__q_h
                          ,"dg " ++ show dg ++ ", call = " ++ show (db_call x)
                          ,"p(h|x) " ++ show p_h__x
                          ,"k  " ++ show k
                          ,"k' " ++ show k'
                          ,"meh " ++ show meh]) $ -} l' gt acc' k' xs

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

{- showCall :: (a -> ShowS) -> VarCall (GL,a) -> ShowS
showCall f vc = shows (vc_refseq vc) . (:) ':' .
                shows (vc_pos vc) . (:) '\t' .
                f (snd $ vc_vars vc) . (++) "\tDP=" .
                shows (vc_depth vc) . (++) ":MQ0=" .
                shows (vc_mapq0 vc) . (++) ":MAPQ=" .
                shows mapq . (:) '\t' .
                show_pl (fst $ vc_vars vc)
  where
    show_pl :: Vector Prob -> ShowS
    show_pl = (++) . intercalate "," . map show . V.toList

    mapq = vc_sum_mapq vc `div` vc_depth vc -}


-- | Error model with dependency parameter.  Since both strands are
-- supposed to still be independent, we feed in only one pile, and
-- later combine both calls.  XXX What's that doing HERE?!

