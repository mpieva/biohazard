{-# LANGUAGE RecordWildCards, OverloadedStrings, FlexibleInstances, BangPatterns #-}
{-# LANGUAGE FlexibleContexts                                                    #-}

-- | Code to read a BAM file in.  We pileup, then sample a base
-- randomly.  We end with the same format we would get if we ran 'vcflc'
-- on an appropriately generated VCF file.

import Bio.Adna                     ( scalarMat )
import Bio.Bam               hiding ( Ns )
import Bio.Bam.Pileup
import Bio.Base              hiding ( Ns )
import Bio.Prelude           hiding ( Ns )
import Data.Bits
import Data.ByteString.Builder      ( hPutBuilder )
import Data.Char
import Data.Foldable                ( toList )
import System.Console.GetOpt
import System.IO
import System.Random                ( randomRIO )
import System.Environment

import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy as L
import qualified Data.Vector.Unboxed as U

-- import Stretch
-- import VcfScan                      ( hashChrom, RawVariant(..) )
-- import Util

main :: IO ()
main = do
    bams <- getArgs
    mergeInputs combineCoordinates bams >=> run $ \hdr ->
            takeWhileE (isValidRefseq . b_rname . unpackBam)            =$
            mapMaybeStream (decompose (repeat (scalarMat 1)))           =$
            pileup                                                      =$
            -- mapStreamM (pick $ chromtab hdr)                            =$
            -- importBam chroms                                            =$
            skipToEof
            -- mapChunksM_ (hPutBuilder hdl . encode_v0 . ($ Done))
  -- where
    -- chromtab  = U.fromList . map hashChrom . map sq_name . toList . meta_refs

-- Huh?  Shouldn't the ref allele be around somewhere?
{-
pick :: U.Vector Int -> Pile -> IO RawVariant
pick chs Pile{ p_snp_pile = (fpile,rpile), .. } = do
    let bases = map (db_call . snd) $ fpile ++ rpile

    i <- randomRIO (0, length bases -1)
    let refbase = db_ref . snd . head $ fpile ++ rpile
        thebase = bases !! i

    return $ if refbase == nucToNucs thebase
        then RawVariant { rv_chrom = chs U.! fromIntegral (unRefseq p_refseq)
                        , rv_pos   = p_pos+1
                        , rv_vars  = B.singleton (showNucleotides refbase)
                        , rv_gt    = 0xFF02 }
        else RawVariant { rv_chrom = chs U.! fromIntegral (unRefseq p_refseq)
                        , rv_pos   = p_pos+1
                        , rv_vars  = showNucleotides refbase `B.cons` (',' `B.cons` B.singleton (showNucleotide thebase))
                        , rv_gt    = 0xFF04 }

importBam :: Monad m => [ L.ByteString ] -> Enumeratee [ RawVariant ] (Stretch -> Stretch) m b
importBam cs0 out0 = getTwo $ nextChrom cs0 out0
  where
    nextChrom [    ] = \out _ _ -> eneeCheckIfDone (\k -> return $ k $ Chunk Break) out
    nextChrom (c:cs) =             generic cs (hashChrom $ B.concat $ L.toChunks c) 1

    generic cs !hs !pos out v1 v2 = eneeCheckIfDone (\k -> generic' cs hs pos k v1 v2) out

    generic'  _ !_     _ k  Nothing        _ = return $ k $ Chunk Break
    generic' cs !hs !pos k (Just var1) mvar2
        | hs /= rv_chrom var1 = nextChrom cs (k $ Chunk Break) (Just var1) mvar2

        -- long gap, creates Ns
        | rv_pos var1 >= pos + 2 = let l  = (rv_pos var1 - pos) `div` 2
                                       k' = k . Chunk $ Ns (fromIntegral l)
                                   in generic cs hs (pos + 2*l) k' (Just var1) mvar2

        -- small gap, have to emit codes
        | rv_pos var1 == pos + 1 = let k' = k . Chunk $ Chrs (NucCode 0) (get_nuc_code var1)
                                   in tryHead >>= generic cs hs (pos+2) k' mvar2

        -- positions must match now
        | rv_pos var1 == pos = case mvar2 of
                -- two variant calls next to each other
                -- if both are reference, we can try and build a stretch
                Just var2 | rv_pos var2 == pos+1 ->
                    if isVar var1 || isVar var2
                      then let k' = k . Chunk $ Chrs (get_nuc_code var1) (get_nuc_code var2)
                           in getTwo $ generic cs hs (pos+2) k'
                      else getTwo $ matches cs hs 1 (pos+2) k

                -- one followed by gap, will become an N
                _ -> let k' = k . Chunk $ Chrs (get_nuc_code var1) (NucCode 0)
                     in tryHead >>= generic cs hs (pos+2) k' mvar2

        | otherwise = error $ "Got variant position " ++ show (rv_pos var1)
                           ++ " when expecting " ++ show pos ++ " or higher."

    -- To extend a stretch of matches, we need
    -- two non-vars at the next two positions
    matches cs !hs !num !pos k (Just var1) (Just var2)
        | rv_pos var1 == pos && rv_pos var2 == pos+1 && not (isVar var1) && not (isVar var2)
            = getTwo $ matches cs hs (num+1) (pos+2) k

    -- anything else, we dump the stretch out and pass the buck
    matches cs hs num pos k v1 v2 = generic cs hs pos (k . Chunk $ Eqs num) v1 v2

    getTwo k = tryHead >>= \a -> tryHead >>= \b -> k a b

    -- *sigh*  Have to turn a numeric genotype into a 'NucCode'.  We
    -- have characters for the variants, and we need to map a pair of
    -- them to a code.
    get_nuc_code RawVariant{..}
        | rv_gt .&. 0xFF00 == 0xFF00        -- haploid call
            = let z = toUpper $ safeIndex "z" rv_vars ( fromIntegral (rv_gt .&. 0x00FE - 2) )
              in NucCode $ maybe (error $ "What's a " ++ shows z "?") fromIntegral
                         $ B.elemIndex z nuc_chars
        | otherwise                         -- diploid
            = let c1 = toUpper $ safeIndex "c1" rv_vars ( fromIntegral (rv_gt            .&. 0x00FE - 2) )
                  c2 = toUpper $ safeIndex "c2" rv_vars ( fromIntegral (rv_gt `shiftR` 8 .&. 0x00FE - 2) )

                  n1 = maybe (error $ "What's a " ++ shows c1 "?") id $ B.elemIndex c1 nuc_chars
                  n2 = maybe (error $ "What's a " ++ shows c2 "?") id $ B.elemIndex c2 nuc_chars
              in NucCode $ two_to_code U.! (n1+5*n2)

    two_to_code = U.fromList [0,1,2,3,4,1,1,5,6,7,2,5,2,8,9,3,6,8,3,10,4,7,9,10,4]

    safeIndex m s i | B.length s > i = B.index s i
                    | otherwise = error $ "Attempted to index " ++ shows i " in " ++ shows s " (" ++ m ++ ")."

    nuc_chars :: B.ByteString
    nuc_chars = "NACGT"

    isVar RawVariant{..} | rv_gt == 0xFF02            = False     -- "0"
                         | rv_gt .&. 0xFCFE == 0x0002 = False     -- "0|.", "0/.", "0|0", "0/0"
                         | otherwise                  = True
-}
