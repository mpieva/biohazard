{-# LANGUAGE BangPatterns #-}
-- Re-wrap alignments to obey the given length of the reference
-- sequence.
--
-- The idea is that a circular reference sequence has been extended
-- artificially to facilitate alignment.  Now the declared length in the
-- header is wrong, and the alignments overhang the end.  Here we split
-- those alignments into two, one for the beginning, one for the end of
-- the sequence, then soft-mask out the inappropriate parts.
--
-- What's the best course of action, operationally?  As usual, we need
-- to decide whether to rely on sorted input and whether to produce
-- sorted output, and how much to copy senselessly.
--
-- In a sane world, this program runs precisely once, after alignment,
-- and output is piped somewhere.  So that's what we do:  input is
-- unsorted, so is output, output is piped (and hence uncompressed).
-- We also fix the header while we're at it.
--
-- We try to fix the map quality for the affected reads as follows:  if
-- a read has map quality 0 (meaning multiple equally good hits), we
-- check the XA field.  If it reports exactly one additional alignment,
-- and it matches the primary alignment when transformed to canonical
-- coordinates, we remove XA and set MAPQ to 37.

import Bio.Bam
import Bio.Base
import Control.Monad                    ( when )
import Paths_biohazard_tools            ( version )
import System.Environment               ( getArgs, getProgName )
import System.Exit                      ( exitFailure )
import System.IO                        ( hPutStr )

import qualified Data.ByteString.Char8  as S
import qualified Data.Map               as M
import qualified Data.Sequence          as Z

usage :: IO a
usage = do pn <- getProgName
           hPutStr stderr $ "Usage: " ++ pn ++ "[chrom:length...]\n\
                \Pipes a BAM file from stdin to stdout and for every 'chrom'\n\
                \mentioned on the command line, wraps alignments to a new \n\
                \target length of 'length'.\n"
           exitFailure

main :: IO ()
main = getArgs >>= \args ->
       when (null args) usage >>= \_ ->
       enumHandle defaultBufSize stdin >=> run $
       joinI $ decodeAnyBam $ \hdr -> do
           add_pg <- liftIO (addPG $ Just version)
           let (ltab, seqs') = parseArgs (meta_refs hdr) args
           joinI $ mapChunks (concatMap (rewrap (M.fromList ltab)))
                 $ pipeRawBamOutput (add_pg hdr { meta_refs = seqs' })

parseArgs :: Refs -> [String] -> ([(Refseq,(Int,S.ByteString))], Refs)
parseArgs = foldl parseArg . (,) []
  where
    parseArg (sqs, h) arg = case break (==':') arg of
        (nm,':':r) -> case reads r of
            [(l,[])] | l > 0 -> case Z.findIndexL ((==) nm . unpackSeqid . sq_name) h of
                Just k  -> case h `Z.index` k of
                    a | sq_length a >= l -> ( (Refseq $ fromIntegral k,(l, sq_name a)):sqs, Z.update k (a { sq_length = l }) h )
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


-- We normalize matches, deletions and soft masks, because these are the
-- operations we generate.  Everything is either already normalized or
-- nobody really cares anyway.
toCigar :: ECig -> Cigar
toCigar = Cigar . go
  where
    go Empty = []
    go (Ins' n ecs) = (Ins,n) : go ecs
    go (Nop' n ecs) = (Nop,n) : go ecs
    go (HMa' n ecs) = (HMa,n) : go ecs
    go (Pad' n ecs) = (Pad,n) : go ecs
    go (SMa' n ecs) = go_sma n ecs
    go (Mat' n ecs) = go_mat n ecs
    go (Rep' _ ecs) = go_mat 1 ecs
    go (Del' _ ecs) = go_del 1 ecs

    go_sma !n (SMa' m ecs) = go_sma (n+m) ecs
    go_sma !n         ecs  = (SMa,n) : go ecs

    go_mat !n (Mat' m ecs) = go_mat (n+m) ecs
    go_mat !n (Rep' _ ecs) = go_mat (n+1) ecs
    go_mat !n         ecs  = (Mat,n) : go ecs

    go_del !n (Del' _ ecs) = go_del (n+1) ecs
    go_del !n         ecs  = (Del,n) : go ecs



setMD :: BamRec -> ECig -> BamRec
setMD b ec = if any interesting md then b { b_exts = M.insert "MD" (Text $ showMd md) (b_exts b) }
                                   else b { b_exts = M.delete "MD" (b_exts b) }
  where
    md = norm $ go ec

    go  Empty       = []
    go (Ins' _ ecs) = go ecs
    go (Nop' _ ecs) = go ecs
    go (SMa' _ ecs) = go ecs
    go (HMa' _ ecs) = go ecs
    go (Pad' _ ecs) = go ecs
    go (Mat' n ecs) = MdNum  n  : go ecs
    go (Rep' x ecs) = MdRep  x  : go ecs
    go (Del' x ecs) = MdDel [x] : go ecs

    norm (MdNum n : MdNum m : mds) = norm $ MdNum (n+m) : mds
    norm (MdDel u : MdDel v : mds) = norm $ MdDel (u++v) : mds
    norm                (op : mds) = op : norm mds
    norm                        [] = []

    interesting (MdRep n) = n /= gap
    interesting (MdDel ns) = all (/= gap) ns
    interesting _ = False


rewrap :: M.Map Refseq (Int,S.ByteString) -> BamRaw -> [BamRaw]
rewrap m br = case M.lookup (br_rname br) m of
    Just (l,nm) | overhangs l    -> do_wrap nm l
                | otherwise      -> do_simple_wrap nm l
    _                            -> [br]
  where
    overhangs l = not (br_isUnmapped br) && br_pos br < l
                  && l < br_pos br + br_aln_length br

    dups_are_fine nm l = br_mapq br == 0 && all_match_XA (br_extAsString "XA" br)
      where
        all_match_XA s = case S.split ';' s of [xa1,xa2] | S.null xa2 -> one_match_XA xa1 ; _ -> False
        one_match_XA s = case S.split ',' s of (sq:pos:_) | sq == nm -> pos_match_XA pos ; _ -> False
        pos_match_XA s = case S.readInt s of Just (p,z) | S.null z -> int_match_XA p ; _ -> False
        int_match_XA p | p >= 0    =  (p-1) `mod` l == br_pos br `mod` l && not (br_isReversed br)
                       | otherwise = (-p-1) `mod` l == br_pos br `mod` l && br_isReversed br

    do_wrap nm l = let !b = decodeBamEntry br ; dp = dups_are_fine nm l in
                   case split_ecig (l - b_pos b) $ toECig (b_cigar b) (maybe [] id $ getMd b) of
                       (left,right) -> [ encodeBamEntry $ b { b_cigar = toCigar  left
                                                            , b_mpos = b_mpos b `mod` l
                                                            , b_mapq = if dp then 37 else b_mapq b
                                                            , b_exts = (if dp then M.delete "XA" else id) (b_exts b) }
                                                          `setMD` left
                                       , encodeBamEntry $ b { b_cigar = toCigar right
                                                            , b_pos = 0
                                                            , b_mpos = b_mpos b `mod` l
                                                            , b_mapq = if dp then 37 else b_mapq b
                                                            , b_exts = (if dp then M.delete "XA" else id) (b_exts b) }
                                                          `setMD` right ]

    do_simple_wrap nm l = [ mutateBamRaw br $ do setPos (br_pos br `mod` l)
                                                 setMpos (br_mpos br `mod` l)
                                                 when (dups_are_fine nm l) (setMapq 37 >> removeExt "XA") ]

-- | Split an 'ECig' into two at some position.  The position is counted
-- in terms of the reference (therefore, deletions count, insertions
-- don't).  The parts that would be skipped if we were splitting lists
-- are replaced by soft masks.
split_ecig :: Int -> ECig -> (ECig, ECig)
split_ecig _ Empty = (Empty,      Empty)
split_ecig 0   ecs = (mask_all ecs, ecs)

split_ecig i (Ins' n ecs) = case split_ecig i ecs of (u,v) -> (Ins' n u, SMa' n v)
split_ecig i (SMa' n ecs) = case split_ecig i ecs of (u,v) -> (SMa' n u, SMa' n v)
split_ecig i (HMa' n ecs) = case split_ecig i ecs of (u,v) -> (HMa' n u, HMa' n v)
split_ecig i (Pad' n ecs) = case split_ecig i ecs of (u,v) -> (Pad' n u,        v)

split_ecig i (Mat' n ecs)
    | i >= n    = case split_ecig (i-n) ecs of (u,v) -> (Mat' n u, SMa' n v)
    | otherwise = (Mat' i $ SMa' (n-1) $ mask_all ecs, SMa' i $ Mat' (n-i) ecs)

split_ecig i (Rep' x ecs) = case split_ecig (i-1) ecs of (u,v) -> (Rep' x u, SMa' 1 v)
split_ecig i (Del' x ecs) = case split_ecig (i-1) ecs of (u,v) -> (Del' x u,        v)

split_ecig i (Nop' n ecs)
    | i >= n    = case split_ecig (i-n) ecs of (u,v) -> (Nop' n u,        v)
    | otherwise = (Nop' i $ mask_all ecs, Nop' (n-i) ecs)

mask_all :: ECig -> ECig
mask_all Empty = Empty
mask_all (Nop' _ ec) =          mask_all ec
mask_all (HMa' _ ec) =          mask_all ec
mask_all (Pad' _ ec) =          mask_all ec
mask_all (Del' _ ec) =          mask_all ec
mask_all (Rep' _ ec) = SMa' 1 $ mask_all ec
mask_all (Mat' n ec) = SMa' n $ mask_all ec
mask_all (Ins' n ec) = SMa' n $ mask_all ec
mask_all (SMa' n ec) = SMa' n $ mask_all ec

