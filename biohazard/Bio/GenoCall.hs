{-# LANGUAGE BangPatterns #-}
module Bio.GenoCall where

import Bio.Adna
import Bio.Bam.Pileup
import Bio.Bam.Rec
import Bio.Base
import Control.Applicative
import Data.Bits ( testBit )
import Data.Foldable hiding ( sum, product )
import Data.List ( tails, intercalate, sortBy )
import Data.Ord

import qualified Data.Set               as Set
import qualified Data.Vector.Unboxed    as V

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

simple_indel_call :: Int -> IndelPile -> (GL, IndelVars)
simple_indel_call ploidy vars = (simple_call ploidy mkpls vars, vars')
  where
    vars' = Set.toList . Set.fromList $ [ (d, map call i) | (_q,(d,i)) <- vars ]
    call (DB a c g t _) = snd $ maximumBy (comparing fst) [(a,nucA), (c,nucC), (g,nucG), (t, nucT)]

    match = zipWith $ \ (DB a c g t q) (N n) -> let p = sum [ if testBit n 0 then a else 0
                                                            , if testBit n 1 then c else 0
                                                            , if testBit n 2 then g else 0
                                                            , if testBit n 3 then t else 0 ]
                                                    p' = fromQual q
                                                in toProb $ p + p' - p * p'

    mkpls (q,(d,i)) = let !q' = qualToProb q
                      in [ if d /= dr || length i /= length ir
                           then q' else q' + product (match i ir) | (dr,ir) <- vars' ]

-- | Naive SNP call; essentially the GATK model.  We create a function
-- that computes a likelihood for a given base, then hand over to simple
-- call.  Since everything is so straight forward, this works even in
-- the face of damage.

simple_snp_call :: Int -> BasePile -> GL
simple_snp_call ploidy vars = simple_call ploidy mkpls vars
  where
    mkpls (q, DB a c g t qq) = [ toProb $ x + pe*(s-x) | x <- [a,c,g,t] ]
      where
        !p1 = fromQual q
        !p2 = fromQual qq
        !pe = p1 + p2 - p1*p2
        !s  = (a+c+g+t) / 4

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

simple_call :: Int -> (a -> [Prob]) -> [a] -> GL
simple_call ploidy pls = foldl1' (V.zipWith (*)) . map step
  where
    foldl1' _ [    ] = V.singleton 1
    foldl1' f (a:as) = foldl' f a as

    norm = toProb (fromIntegral ploidy) `pow` (-1)

    -- "For biallelic sites the ordering is: AA,AB,BB; for triallelic
    -- sites the ordering is: AA,AB,BB,AC,BC,CC, etc."
    --
    -- To get the order right, we reverse the list, recurse, and reverse
    -- the result again.
    step = V.fromList . map (* norm) . reverse . mk_pls ploidy . reverse . pls

    -- Meh.  Pointless, but happens to be the unit.
    mk_pls 0  _ = return 0

    -- Okay, we sample ONE allele.  Likelihood of the data is simply the
    -- GL value that was passed to us.
    mk_pls 1 ls = ls

    -- We extend the genotype and sample another allele.
    mk_pls n ls = do ls'@(hd:_) <- tails ls
                     (+) hd <$> mk_pls (n-1) ls'


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
    show_pl :: V.Vector Prob -> ShowS
    show_pl = (++) . intercalate "," . map show . V.toList

    mapq = vc_sum_mapq vc `div` vc_depth vc -}


-- | Error model with dependency parameter.  Since both strands are
-- supposed to still be independent, we feed in only one pile, and
-- later combine both calls.

maq_snp_call :: Int -> Double -> BasePile -> GL
maq_snp_call ploidy theta bases = undefined
  where
    bases' = sortBy (\(DB _ _ _ _ q1) (DB _ _ _ _ q2) -> compare q2 q1)
             [ DB a c g t (min q mq) | (mq, DB a c g t q) <- bases ]

