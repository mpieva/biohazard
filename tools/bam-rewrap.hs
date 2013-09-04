-- Re-wrap alignments to obey the given length of the reference
-- sequence.
--
-- The idea is that a circular reference sequence has been extended
-- artificially to facilitate alignment.  Now the declared length in the
-- header is wrong, and the alignments overhang the end.  Here we split
-- those alignments into two, one for the beginning, one for the end of
-- the sequence, then mask out the inappropriate parts.
--
-- What's the best course of action, operationally?  As usual, we need
-- to decide whether to rely on sorted input and whether to produce
-- sorted output, and how much to copy senselessly.
--
-- In a sane world, this program runs precisely once, after alignment,
-- and output is piped somewhere.  So that's what we do:  input is
-- unsorted, so is output, output is piped (and hence uncompressed).
-- We also fix the header while we're at it.

import Bio.Bam
import Bio.Base
import Paths_biohazard_tools            ( version )
import System.Environment               ( getArgs )

import qualified Data.Map               as M
import qualified Data.Sequence          as Z

main :: IO ()
main = enumHandle defaultBufSize stdin >=> run $
       joinI $ decodeAnyBam $ \hdr -> do
           add_pg       <- liftIO (addPG $ Just version)
           (ltab, seqs') <- parseArgs (meta_refs hdr) `fmap` liftIO getArgs
           joinI $ mapChunks (concatMap (rewrap (M.fromList ltab)))
                 $ pipeRawBamOutput (add_pg hdr { meta_refs = seqs' })

parseArgs :: Refs -> [String] -> ([(Refseq,Int)], Refs)
parseArgs = foldl parseArg . (,) []
  where
    parseArg (sqs, h) arg = case break (==':') arg of
        (nm,':':r) -> case reads r of
            [(l,[])] | l > 0 -> case Z.findIndexL ((==) nm . unpackSeqid . sq_name) h of
                Just k  -> case h `Z.index` k of
                    a | sq_length a >= l -> ( (Refseq $ fromIntegral k,l):sqs, Z.update k (a { sq_length = l }) h )
                      | otherwise -> error $ "cannot wrap " ++ show nm ++ " to " ++ show l
                                          ++ ", which is more than the original " ++ show (sq_length a)
                Nothing -> error $ "target sequence " ++ show nm ++ " not found"
            _ -> error $ "couldn't parse length " ++ show r ++ " for " ++ show nm
        _ -> error $ "couldn't parse argument " ++ show arg

-- | Argh, this business with the CIGAR operations is a mess, it gets
-- worse when combined with MD.  Okay, we will support CIGAR (no "=" and
-- "X" operations) and MD.  If we have MD on input, we generate it on
-- output, too.  And in between, we break everything into /very small/
-- operations.

data ECig = Empty
          | Mat' Int ECig
          | Rep' Nucleotide ECig
          | Ins' Int ECig
          | Del' Nucleotide ECig
          | Nop' Int ECig
          | SMa' Int ECig
          | HMa' Int ECig
          | Pad' Int ECig

toECig :: Cigar -> [MdOp] -> ECig
toECig (Cigar cig) md = go cig md
  where
    go [       ]             _  = Empty
    go        cs (MdNum  0:mds) = go cs mds
    go        cs (MdDel []:mds) = go cs mds
    go ((_,0):cs)          mds  = go cs mds

    go ((Mat,n):cs) [           ]      = Mat'   n  $ go            cs              [ ]
    go ((Mat,n):cs) (MdRep x:mds)      = Rep'   x  $ go ((Mat,n-1):cs)             mds
    go ((Mat,n):cs) (MdDel z:mds)      = Mat'   n  $ go            cs     (MdDel z:mds)
    go ((Mat,n):cs) (MdNum m:mds)
       | n < m                         = Mat'   n  $ go            cs (MdNum (m-n):mds)
       | n > m                         = Mat'   m  $ go ((Mat,n-m):cs)             mds
       | otherwise                     = Mat'   n  $ go            cs              mds

    go ((Ins,n):cs)               mds  = Ins'   n  $ go            cs              mds
    go ((Del,n):cs) (MdDel (x:xs):mds) = Del'   x  $ go ((Del,n-1):cs)   (MdDel xs:mds)
    go ((Del,n):cs)               mds  = Del' nucN $ go ((Del,n-1):cs)             mds

    go ((Nop,n):cs) mds = Nop' n $ go cs mds
    go ((SMa,n):cs) mds = SMa' n $ go cs mds
    go ((HMa,n):cs) mds = HMa' n $ go cs mds
    go ((Pad,n):cs) mds = Pad' n $ go cs mds

toCigar :: ECig -> Cigar
toCigar = undefined -- !!!

toMD :: ECig -> [MdOp]
toMD = undefined -- !!!

setMD :: BamRec -> [MdOp] -> BamRec
setMD b [] = b

rewrap :: M.Map Refseq Int -> BamRaw -> [BamRaw]
rewrap m br = case M.lookup (br_rname br) m of
    Just l | not (br_isUnmapped br) && l < br_pos br + br_aln_length br -> do_wrap l
    _                                                                   -> [br]
  where
    b = decodeBamEntry br
    do_wrap l = case split_ecig (l - b_pos b) $ toECig (b_cigar b) (maybe [] id $ getMd b) of
                    (left,right) -> [ encodeBamEntry $ b { b_cigar = toCigar  left } `setMD` toMD  left
                                    , encodeBamEntry $ b { b_cigar = toCigar right } `setMD` toMD right ]

    split_ecig n _ = (Empty,Empty) -- !!!



