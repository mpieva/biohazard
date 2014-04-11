{-# LANGUAGE OverloadedStrings, BangPatterns, RecordWildCards #-}
{-# OPTIONS_GHC -Wall #-}

-- Cobble up a mitochondrion, or something similar.
--
-- The goal is to reconstruct a mitochondrion (or similar small, haploid
-- locus) from a set of sequencing reads and a reference sequence.  The
-- idea is to first select reads using some sort of filtering strategy,
-- simply for speed reasons.  They are then aligned to the reference
-- using banded Smith-Waterman algorithm, and a more likely reference is
-- called.  This is repeated till it converges.  A bad implementation of
-- the idea was called MIA.

import Align
import SimpleSeed

import Bio.Base
import Bio.Bam
import Bio.Iteratee
import Control.Applicative
import Control.Monad
import Data.Char
import Data.Monoid
import Numeric
import System.Console.GetOpt
import System.Environment
import System.Exit
import System.IO

import qualified Data.ByteString            as B
import qualified Data.ByteString.Char8      as S
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.Foldable              as F
import qualified Data.Iteratee              as I
import qualified Data.Map                   as M
import qualified Data.Sequence              as Z
import qualified Data.Vector.Unboxed        as U


-- Read a FastA file, drop the names, yield the sequences.
readFasta :: L.ByteString -> [( S.ByteString, [Either Nucleotide Nucleotide] )]
readFasta = go . dropWhile (not . isHeader) . L.lines
  where
    isHeader s = not (L.null s) && L.head s == '>'

    go [     ] = []
    go (hd:ls) = case break isHeader ls of
                (body, rest) -> let ns = map toNuc . concat $ map L.unpack body
                                    nm = S.concat . L.toChunks . L.tail . head $ L.words hd
                                in (nm,ns) : if null rest then [] else go rest

    toNuc x | isUpper x = Right $ toNucleotide x
            | otherwise = Left  $ toNucleotide (toUpper x)


-- | A query record.  We construct these after the seeding phase and
-- keep the bare minimum:  name, sequence/quality, seed region, flags
-- (currently only the strand).  Just enough to write a valig BAM file.

data QueryRec = QR { qr_name :: {-# UNPACK #-} !Seqid           -- from BAM
                   , qr_seq  :: {-# UNPACK #-} !QuerySeq        -- sequence and quality
                   , qr_pos  :: {-# UNPACK #-} !RefPosn         -- start position of band
                   , qr_band :: {-# UNPACK #-} !Bandwidth }     -- bandwidth (negative to indicate reversed sequence_
  deriving Show

data Conf = Conf {
    conf_references :: [FilePath] -> [FilePath],
    conf_aln_outputs :: Maybe (Int -> FilePath),
    conf_cal_outputs :: Maybe (Int -> FilePath) }

iniconf :: Conf
iniconf = Conf id Nothing Nothing

options :: [ OptDescr (Conf -> IO Conf) ]
options = [
    Option "r" ["reference"]  (ReqArg add_ref "FILE") "Read references from FILE",
    Option "a" ["align-out"]  (ReqArg set_aln_out "PAT") "Write intermediate alignments to PAT",
    Option "c" ["called-out"] (ReqArg set_cal_out "PAT") "Write called references to PAT" ]
  where
    add_ref     f c = return $ c { conf_references = conf_references c . (:) f }
    set_aln_out p c = return $ c { conf_aln_outputs = Just (splice_pat p) }
    set_cal_out p c = return $ c { conf_cal_outputs = Just (splice_pat p) }

    splice_pat [] _ = []
    splice_pat ('%':'%':s) x = '%' : splice_pat s x
    splice_pat ('%':'d':s) x = shows x $ splice_pat s x
    splice_pat (c:s) x = c : splice_pat s x

main :: IO ()
main = do
    (opts, files, errors) <- getOpt Permute options <$> getArgs
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Conf{..} <- foldl (>>=) (return iniconf) opts

    inputs@((refname,reference):_) <- concatMap readFasta <$> mapM L.readFile (conf_references [])

    let !sm = create_seed_maps (map (map (either id id) . snd) inputs)
        !rs = prep_reference reference

    let bamhdr = mempty { meta_hdr = BamHeader (1,4) Unsorted []
                        , meta_refs = Z.singleton $ BamSQ refname (length reference) [] }


    -- uhh.. termination condition?
    let round n k     = do let bamout = case conf_aln_outputs of
                                            Nothing -> skipToEof
                                            Just nf -> write_iter_bam (nf n) bamhdr
                           (newref, queries) <- k bamout
                           case conf_cal_outputs of Nothing -> return ()
                                                    Just nf -> write_ref_fasta (nf n) n newref
                           putStrLn $ "Round " ++ shows n ": Kept " ++ shows (length queries) " queries."
                           round (n+1) (\out -> enumPure1Chunk queries >=> run $ roundN undefined out)

    round 1 (\out -> concatInputs files >=> run $ \_ -> round1 sm rs out)


    -- print queries
    return ()
  where

{- XXX
add_to_refseq' :: MonadIO m => NewRefSeq -> QuerySeq -> AlignResult -> m NewRefSeq
add_to_refseq' nrs@(NRS v) qs res = do
    liftIO $ print (viterbi_position res, viterbi_score res, viterbi_backtrace res)
    liftIO $ print $ take 200 $ F.foldr foo [] v
    liftIO $ print $ take 200 $ ref_to_ascii $ finalize_ref_seq nrs
    let nrs' = add_to_refseq nrs qs res
    liftIO $ print {- $ take 200 -} $ filter (not . U.all (==0)) $ F.foldr foo [] v
    liftIO $ print $ take 200 $ ref_to_ascii $ finalize_ref_seq nrs'
    return nrs'
  where
    foo (NC is b) l = is : b : l


test :: AlignResult
test = align 5 rs qs (RP 2) (BW 3)
  where
    qs = prep_query_fwd $ encodeBamEntry $ nullBamRec { b_seq = U.fromList $ map toNucleotide "ACGT", b_qual = B.pack [20,21,22,23] }
    -- qs = prep_query_fwd $ encodeBamEntry $ nullBamRec { b_seq = U.fromList $ map toNucleotide "AC", b_qual = B.pack [20,21] }
    rs = prep_reference $ map Right . map toNucleotide $ "AAAACCGTTTT"
-}

-- General plan:  In the first round, we read, seed, align, call the new
-- working sequence, and write a BAM file.  Then write the new working
-- sequence out.  In subsequent rounds, the seeding is skipped and the
-- sequences come from memory.
--
-- XXX the bandwidth is too low, definitely in round 1, probably in
-- subsequent rounds.

round1 :: MonadIO m
       => SeedMap -> RefSeq
       -> Iteratee [(QueryRec, AlignResult)] m ()        -- BAM output
       -> Iteratee [BamRaw] m         -- queries in
            (RefSeq, [QueryRec])  -- new reference & queries out
round1 sm rs out = convStream (headStream >>= seed) =$ roundN rs out
  where
    seed br | low_qual br = return []
    seed br = case do_seed (refseq_len rs) sm br of
        Nothing    -> return []
        Just (a,b) -> let bw = b - a - br_l_seq br
                      in return $ if a >= 0
                         then [ QR (br_qname br) (prep_query_fwd br) (RP   a ) (BW   bw ) ]
                         else [ QR (br_qname br) (prep_query_rev br) (RP (-b)) (BW (-bw)) ]

    low_qual br = 2 * l1 > l2 where
        l2 = br_l_seq br
        l1 = F.foldl' step 0 [0 .. l2-1]
        step a i = if br_qual_at br i > Q 10 then a+1 else a

roundN :: Monad m
       => RefSeq
       -> Iteratee [(QueryRec, AlignResult)] m ()        -- BAM output
       -> Iteratee [QueryRec] m         -- queries in
            (RefSeq, [QueryRec])  -- new reference & queries out
roundN rs out = do
    ((), (rs', xtab), qry') <- mapStream aln =$ filterStream good =$
                               I.zip3 out mkref collect
    return (rs', reverse $ map (xlate xtab) qry')

  where
    gap_cost = 50           -- Hmm, better suggestions?
    pad = 8

    aln qr@QR{..} = let res = align gap_cost rs qr_seq qr_pos qr_band
                    in ( new_coords qr res, res )
    good (_, res) = viterbi_score res < 0

    mkref = finalize_ref_seq `liftM` foldStream step (new_ref_seq rs)
    step nrs (qr, res) = add_to_refseq nrs (qr_seq qr) res

    collect :: Monad m => Iteratee [(QueryRec, AlignResult)] m [(Int,Int,QueryRec)]
    collect = foldStream (\l (!qr,!ar) ->
                -- get alignment ends from ar, add some buffer
                let !left  = viterbi_position ar - 8
                    !right = viterbi_position ar + 8 + cigarToAlnLen (viterbi_backtrace ar)
                in (left,right,qr) : l) []

    xlate :: XTab -> (Int, Int, QueryRec) -> QueryRec
    xlate tab (l,r,qr) =
        let !left  = Z.index tab l
            !right = Z.index tab r
        in qr { qr_pos = RP left, qr_band = BW $ right - left }

    new_coords qr rs = qr { qr_pos  = RP $ viterbi_position rs - pad
                          , qr_band = BW $ (if reversed (qr_band qr) then negate else id) $
                                      2*pad + max_bandwidth (viterbi_backtrace rs) }

    reversed (BW x) = x < 0

    max_bandwidth = (+1) . (*2) . maximum . map abs . scanl plus 0 . unCigar
    plus a (Mat,_) = a
    plus a (Ins,n) = a+n
    plus a (Del,n) = a-n



-- Outline for further rounds:  We keep the same queries, we use the new
-- reference called in the previous round.  Output is channelled to
-- different files.  However, we need to translate coordinates to keep
-- the alignment windows in the correct places.  This should actually
-- come from the calling of the new reference.  Note that coordinate
-- translation may actually change the bandwidth.  Also we have to
-- compute a sensible bandwidth from the alignment.

write_iter_bam :: FilePath -> BamMeta -> Iteratee [(QueryRec, AlignResult)] IO ()
write_iter_bam fp hdr = mapStream conv =$ writeBamFile fp hdr
  where
    conv (QR{..}, AlignResult{..}) = BamRec
            { b_qname           = qname
            , b_flag            = if reversed qr_band then flagReversed else 0
            , b_rname           = Refseq 0
            , b_pos             = viterbi_position
            , b_mapq            = 255
            , b_cigar           = viterbi_backtrace
            , b_mrnm            = invalidRefseq
            , b_mpos            = 0
            , b_isize           = 0
            , b_seq             = qseqToBamSeq qr_seq
            , b_qual            = qseqToBamQual qr_seq
            , b_virtual_offset  = 0
            , b_exts            = M.empty }
      where
        qname = qr_name `S.append` S.pack ("  " ++ showFFloat (Just 1) viterbi_score [])
        reversed (BW x) = x < 0

-- Call sequence and write to file.  We call a base only if all bases
-- together are more likely than a gap.  We call a weak base if the gap
-- has a probality of more than 30%.  The called base is the most
-- probable one.
write_ref_fasta :: FilePath -> Int -> RefSeq -> IO ()
write_ref_fasta fp num rs = writeFile fp $ unlines $
    (">genotype_call-" ++ show num) : chunk 70 (ref_to_ascii rs)
  where
    chunk n s = case splitAt n s of _ | null s -> [] ; (l,r) -> l : chunk n r

ref_to_ascii :: RefSeq -> String
ref_to_ascii (RS v) = [ base | i <- [0, 5 .. U.length v - 5]
                             , let pgap = v U.! (i+4)
                             , pgap <= 3
                             , let letters = if pgap >= 1 then "acgtn" else "ACGTN"
                             , let index = U.maxIndex (U.slice i 4 v)
                             , let good = U.maximum (U.slice i 4 v) > 3
                             , let base = S.index letters $ if good then index else 4 ]

