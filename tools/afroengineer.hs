{-# LANGUAGE OverloadedStrings, BangPatterns, RecordWildCards, RankNTypes #-}
{-# OPTIONS_GHC -Wall #-}

-- Cobble up a mitochondrion, or something similar.  This is not an
-- assembly, but something that could serve in stead of one :)
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
import Control.Applicative
import Control.Monad
import Data.Bits
import Data.Char
import Data.List ( isSuffixOf )
import Data.Monoid
import Numeric
import Prelude hiding ( round )
import System.Console.GetOpt
import System.Directory ( doesFileExist )
import System.Environment
import System.Exit
import System.IO

import qualified Bio.Iteratee.ZLib          as ZLib
import qualified Data.ByteString.Char8      as S
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.Foldable              as F
import qualified Data.Iteratee              as I
import qualified Data.Sequence              as Z
import qualified Data.Vector.Unboxed        as U

import Debug.Trace

-- Read a FastA file, drop the names, yield the sequences.
readFasta :: L.ByteString -> [( S.ByteString, [Either Nucleotides Nucleotides] )]
readFasta = go . dropWhile (not . isHeader) . L.lines
  where
    isHeader s = not (L.null s) && L.head s == '>'

    go [     ] = []
    go (hd:ls) = case break isHeader ls of
                (body, rest) -> let ns = map toNuc . concat $ map L.unpack body
                                    nm = S.concat . L.toChunks . L.tail . head $ L.words hd
                                in (nm,ns) : if null rest then [] else go rest

    toNuc x | isUpper x = Right $ toNucleotides x
            | otherwise = Left  $ toNucleotides (toUpper x)


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
                           round (n+1) (\out -> enumPure1Chunk queries >=> run $ roundN newref out)

    round 1 (\out -> foldr ((>=>) . readFreakingInput) run files $ round1 sm rs out)

    -- print queries
    -- return ()


-- General plan:  In the first round, we read, seed, align, call the new
-- working sequence, and write a BAM file.  Then write the new working
-- sequence out.  In subsequent rounds, the seeding is skipped and the
-- sequences come from memory.
--
-- XXX the bandwidth is too low, definitely in round 1, probably in
-- subsequent rounds.

round1 :: MonadIO m
       => SeedMap -> RefSeq
       -> Iteratee [(QueryRec, AlignResult)] m ()       -- BAM output
       -> Iteratee [BamRaw] m                           -- queries in
            (RefSeq, [QueryRec])                        -- new reference & queries out
round1 sm rs out = convStream (headStream >>= seed) =$ roundN rs out
  where
    seed br = case do_seed (refseq_len rs) sm br of
        _ | low_qual br        -> return []
        Nothing                -> return []
        Just (a,b) | a >= 0    -> return [ QR (br_qname br) (prep_query_fwd br) (RP   a ) (BW   bw ) ]
                   | otherwise -> return [ QR (br_qname br) (prep_query_rev br) (RP (-b)) (BW (-bw)) ]
            where bw = b - a - br_l_seq br

    low_qual br = 2 * l1 < l2 where
        l2 = br_l_seq br
        l1 = F.foldl' step 0 [0 .. l2-1]
        step a i = if br_qual_at br i > Q 10 then a+1 else a

roundN :: Monad m
       => RefSeq
       -> Iteratee [(QueryRec, AlignResult)] m ()       -- BAM output
       -> Iteratee [QueryRec] m                         -- queries in
            (RefSeq, [QueryRec])                        -- new reference & queries out
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
                -- XXX does this yield invalid coordinates?
                let !left  = viterbi_position ar - 8
                    !right = viterbi_position ar + 8 + cigarToAlnLen (viterbi_backtrace ar)
                in (left,right,qr) : l) []

    xlate :: XTab -> (Int, Int, QueryRec) -> QueryRec
    xlate tab (l,r,qr)
        | r <= l = error "confused reft and light"
        | left < 0 || right < 0 = error "too far left"
        | right' < left = error "flipped over"
        | otherwise = qr { qr_pos = RP left, qr_band = BW $ right' - left }
      where
        lk x | x < 0            = Z.index tab (x + Z.length tab - 1)
             | x < Z.length tab = Z.index tab  x
             | otherwise        = Z.index tab (x - Z.length tab + 1)

        left = lk l ; right = lk r ; _ Z.:> newlen = Z.viewr tab
        right' = if left < right then right else right + newlen

    new_coords qr rr = qr { qr_pos  = RP $ viterbi_position rr - pad
                          , qr_band = BW $ (if reversed (qr_band qr) then negate else id) $
                                      2*pad + max_bandwidth (viterbi_backtrace rr) }

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
            , b_mapq            = Q 255
            , b_cigar           = viterbi_backtrace
            , b_mrnm            = invalidRefseq
            , b_mpos            = 0
            , b_isize           = 0
            , b_seq             = qseqToBamSeq qr_seq
            , b_qual            = qseqToBamQual qr_seq
            , b_virtual_offset  = 0
            , b_exts            = [] }
      where
        qname = qr_name `S.append` S.pack ("  " ++ showFFloat (Just 1) viterbi_score [])
        reversed (BW x) = x < 0

-- | Calls sequence and writes to file.  We call a base only if the gap
-- has a probability lower than 50%.  We call a weak base if the gap has
-- a probality of more than 25%.  If the most likely base is at least
-- twice as likely as the second most likely one, we call it.  Else we
-- call an N or n.
write_ref_fasta :: FilePath -> Int -> RefSeq -> IO ()
write_ref_fasta fp num rs = writeFile fp $ unlines $
    (">genotype_call-" ++ show num) : chunk 70 (ref_to_ascii rs)
  where
    chunk n s = case splitAt n s of _ | null s -> [] ; (l,r) -> l : chunk n r

ref_to_ascii :: RefSeq -> String
ref_to_ascii (RS v) = [ base | i <- [0, 5 .. U.length v - 5]
                             , let pgap = indexV "ref_to_ascii/pgap" v (i+4)
                             , pgap > 3
                             , let letters = if pgap <= 6 then "acgtn" else "ACGTN"
                             , let (index, p1, p2) = minmin i 4
                             , let good = p2 - p1 >= 3 -- probably nonsense
                             , let base = S.index letters $ if good then index else  trace (show (U.slice i 5 v)) 4 ]
  where
    minmin i0 l = U.ifoldl' step (l, 255, 255) $ U.slice i0 l v
    step (!i, !m, !n) j x | x <= m    = (j, x, m)
                          | x <= n    = (i, m, x)
                          | otherwise = (i, m, n)



readFreakingInput :: (MonadIO m, MonadMask m) => FilePath -> Enumerator [BamRaw] m b
readFreakingInput fp k | ".bam" `isSuffixOf` fp = do liftIO (hPutStrLn stderr $ "Reading BAM from " ++ fp)
                                                     decodeAnyBamFile fp $ const k
                       | otherwise              = maybe_read_two fp unzipFastq k

check_r2 :: FilePath -> IO (Maybe FilePath)
check_r2 = go [] . reverse
  where
    go acc ('1':'r':fp) = do let fp' = reverse fp ++ 'r' : '2' : acc
                             e <- doesFileExist fp'
                             return $ if e then Just fp' else Nothing
    go acc (c:fp) = go (c:acc) fp
    go  _  [    ] = return Nothing

maybe_read_two :: (MonadIO m, MonadMask m)
    => FilePath
    -> (forall m1 b . (MonadIO m1, MonadMask m1) => Enumeratee S.ByteString [BamRec] m1 b)
    -> Enumerator [BamRaw] m a
maybe_read_two fp e1 = (\k -> liftIO (check_r2 fp) >>= maybe (rd1 k) (rd2 k)) $= mapStream encodeBamEntry
  where
    rd1 k     = do liftIO (hPutStrLn stderr $ "Reading FastQ from " ++ fp)
                   enumFile defaultBufSize fp  $= e1 $ k
    rd2 k fp' = do liftIO (hPutStrLn stderr $ "Reading FastQ from " ++ fp ++ " and " ++ fp')
                   mergeEnums (enumFile defaultBufSize fp  $= e1)
                              (enumFile defaultBufSize fp' $= e1)
                              (convStream unite_pairs) k

-- No, we don't need to 'removeWarts'.  This input is, of course, a special case.  :-(
unzipFastq :: (MonadIO m, MonadMask m) => Enumeratee S.ByteString [BamRec] m b
unzipFastq = ZLib.enumInflateAny ><> parseFastq

unite_pairs :: Monad m => Iteratee [BamRec] (Iteratee [BamRec] m) [BamRec]
unite_pairs = do a <- lift headStream
                 b <- headStream
                 return [ a { b_flag = b_flag a .|. flagFirstMate }
                        , b { b_flag = b_flag b .|. flagSecondMate } ]

