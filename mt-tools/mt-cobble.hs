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

    let round1out = case conf_aln_outputs of
                        Nothing -> skipToEof
                        Just nf -> write_iter_bam (nf 1) bamhdr

    ((), newref, queries) <- concatInputs files >=> run $ \_ ->
                                joinI $ round1 sm rs $ I.zip3 round1out (mknewref rs) collect

    case conf_cal_outputs of Nothing -> return ()
                             Just nf -> write_ref_fasta (nf 1) newref

    putStrLn $ "Kept " ++ shows (length queries) " queries."
    return ()
  where
    collect = foldStream (flip (:)) []
    mknewref rs = foldStream (\nrs (qr, res) -> add_to_refseq nrs (qr_seq qr) res)
                             (new_ref_seq rs) >>= return . finalize_ref_seq

test :: AlignResult
test = align 5 rs qs (RP 2) (BW 3)
  where
    qs = prep_query_fwd $ encodeBamEntry $ nullBamRec { b_seq = U.fromList $ map toNucleotide "ACGT", b_qual = B.pack [20,21,22,23] }
    -- qs = prep_query_fwd $ encodeBamEntry $ nullBamRec { b_seq = U.fromList $ map toNucleotide "AC", b_qual = B.pack [20,21] }
    rs = prep_reference $ map Right . map toNucleotide $ "AAAACCGTTTT"

-- General plan:  In the first round, we read, seed, align, call the new
-- working sequence, and write a BAM file.  Then write the new working
-- sequence out.  In subsequent rounds, the seeding is skipped and the
-- sequences come from memory.

round1 :: MonadIO m
       => SeedMap -> RefSeq
       -> Enumeratee [BamRaw] [(QueryRec, AlignResult)] m a
round1 sm rs = convStream (headStream >>= go)
  where
    -- Hmm, better suggestions?
    gap_cost = 50

    go br = case do_seed (refseq_len rs) sm br of
        Nothing    -> return []
        Just (a,b) -> let bw = b - a - br_l_seq br
                      in return $ if a >= 0 then aln (prep_query_fwd br) (RP   a ) (BW bw)
                                            else aln (prep_query_rev br) (RP (-b)) (BW (-bw))
      where
        aln qs (RP x) (BW y) =
            let res = align gap_cost rs qs (RP x) (BW (abs y))
            in if viterbi_score res >= 0 then [] else
               [( QR (br_qname br) qs (RP x) (BW y), res )]

            -- let AlignResult{..} =
            --

write_iter_bam :: MonadCatchIO m => FilePath -> BamMeta -> Iteratee [(QueryRec, AlignResult)] m ()
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
write_ref_fasta :: FilePath -> RefSeq -> IO ()
write_ref_fasta fp (RS v) = writeFile fp $ unlines $
    ">genotype_call" : chunk 60 cseq
  where
    cseq = [ base | i <- [0, 5 .. U.length v - 5]
                  , let pgap = v U.! (i+4)
                  , pgap <= 3
                  , let letters = if pgap >= 1 then "acgt" else "ACGT"
                  , let base = letters `S.index` U.maxIndex (U.slice i 4 v) ]

    chunk n s = case splitAt n s of _ | null s -> [] ; (l,r) -> l : chunk n r
