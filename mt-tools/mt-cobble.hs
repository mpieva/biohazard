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
import System.Console.GetOpt
import System.Environment
import System.Exit
import System.IO

import qualified Data.ByteString.Char8      as S
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.Map                   as M
import qualified Data.Sequence              as Z


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

data QueryRec = QR { qr_name :: {-# UNPACK #-} !Seqid
                   , qr_seq  :: {-# UNPACK #-} !QuerySeq
                   , qr_pos  :: {-# UNPACK #-} !RefPosn
                   , qr_band :: {-# UNPACK #-} !Bandwidth
                   , qr_flag :: {-# UNPACK #-} !Int }

data Conf = Conf {
    conf_references :: [FilePath] -> [FilePath],
    conf_aln_outputs :: Maybe (Int -> FilePath) }

iniconf :: Conf
iniconf = Conf id Nothing

options :: [ OptDescr (Conf -> IO Conf) ]
options = [
    Option "r" ["reference"] (ReqArg add_ref "FILE") "Read references from FILE",
    Option "a" ["align-out"] (ReqArg set_aln_out "PAT") "Write intermediate alignments to PAT" ]
  where
    add_ref     f c = return $ c { conf_references = conf_references c . (:) f }
    set_aln_out p c = return $ c { conf_aln_outputs = Just (splice_pat p) }

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
                        Just nf -> writeBamFile (nf 1) bamhdr

    concatDefaultInputs >=> run >=> run $ \_ -> round1 sm rs round1out



-- General plan:  In the first round, we read, seed, align, call the new
-- working sequence, and write a BAM file.  Then write the new working
-- sequence out.  In subsequent rounds, the seeding is skipped and the
-- sequences come from memory.

round1 :: MonadIO m
       => SeedMap -> RefSeq -> Iteratee [BamRec] m ()
       -> Iteratee [BamRaw] m ( Iteratee [BamRec] m () )
round1 sm rs = convStream (headStream >>= go)
  where
    -- Hmm, better suggestions?
    gap_cost = 50

    go br = case do_seed (refseq_len rs) sm br of
        Nothing    -> return []
        Just (a,b) -> let bw = BW $ b - a - br_l_seq br
                      in if a >= 0 then aln (prep_query_fwd br) (RP   a ) bw
                                   else aln (prep_query_rev br) (RP (-b)) bw
      where
        aln  qs (RP x) (BW y) = do
            -- liftIO $ do putStrLn $ "Rgn " ++ show x ++ ".." ++ show (x+br_l_seq br) ++ "x" ++ show y
                        -- hFlush stdout
            -- liftIO $ print qs
            let memo = viterbi_forward gap_cost rs qs (RP x) (BW y)
            let mm = get_memo_max (BW y) memo
            -- liftIO $ print mm

            return [ (decodeBamEntry br) { b_rname = Refseq 0
                                         -- , b_flag :: !Int
                                         , b_pos = x
                                         , b_mapq = 255
                                         , b_cigar = Cigar []
                                         , b_mrnm = invalidRefseq
                                         , b_mpos = 0
                                         , b_isize = 0
                                         , b_virtual_offset = 0
                                         , b_exts = M.singleton "AS" (Int mm)
                                         } ]
                -- b_seq :: !Bio.Base.Sequence,
                -- b_qual :: !Bio.Bam.Rec.ByteString,
                -- b_exts :: Extensions,

