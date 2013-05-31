{-# LANGUAGE RecordWildCards, BangPatterns #-}
import Bio.File.Bam
import Bio.File.Bam.Rmdup
import Bio.Iteratee
import Bio.Util ( showNum, showOOM, estimateComplexity )
import Control.Monad
import Data.Bits
import Data.List ( intercalate )
import Data.Maybe
import Data.Monoid ( mempty )
import Data.Word ( Word8 )
import Numeric ( showFFloat )
import Paths_biohazard ( version )
import System.Console.GetOpt
import System.Environment ( getArgs, getProgName )
import System.Exit
import System.IO

import qualified Data.ByteString    as S
import qualified Data.Map           as M
import qualified Data.Iteratee      as I

data Conf = Conf {
    output :: BamMeta -> Iteratee [BamRaw] IO (),
    strand_preserved :: Bool,
    collapse :: Bool -> Collapse,
    keep_all :: Bool,
    keep_unaligned :: Bool,
    filter_enee :: BamRaw -> Maybe BamRaw,
    min_len :: Int,
    get_label :: M.Map Seqid Seqid -> BamRaw -> Seqid,
    debug :: String -> IO (),
    which :: Which }

-- | Which reference sequences to scan
data Which = All | Some Refseq Refseq | Unaln deriving Show

defaults :: Conf
defaults = Conf { output = pipeRawBamOutput
                , strand_preserved = True
                , collapse = cons_collapse' 60
                , keep_all = False
                , keep_unaligned = False
                , filter_enee = is_aligned
                , min_len = 0
                , get_label = get_library
                , debug = \_ -> return ()
                , which = All }

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option  "o" ["output"]         (ReqArg set_output "FILE") "Write to FILE (default: stdout)",
    Option  "R" ["refseq"]         (ReqArg set_range "RANGE") "Read only range of reference sequences",
    Option  "p" ["improper-pairs"] (NoArg  set_improper)      "Include improper pairs",
    Option  "u" ["unaligned"]      (NoArg  set_unaligned)     "Included unaligned reads and pairs",
    Option  "1" ["single-read"]    (NoArg  set_single)        "Pretend there is no second mate",
    Option  "c" ["cheap"]          (NoArg  set_cheap)         "Cheap computation: skip the consensus calling",
    Option  "C" ["count-only"]     (NoArg  set_count_only)    "Count duplicates, don't bother with output",
    Option  "k" ["keep","mark-only"](NoArg set_keep)          "Mark duplicates, but include them in output",
    Option  "Q" ["max-qual"]       (ReqArg set_qual "QUAL")   "Set maximum quality after consensus call to QUAL",
    Option  "l" ["min-length"]     (ReqArg set_len "LEN")     "Discard reads shorter than LEN",
    Option  "s" ["no-strand"]      (NoArg  set_no_strand)     "Strand of alignments is uninformative",
    Option  "r" ["ignore-rg"]      (NoArg  set_no_rg)         "Ignore read groups when looking for duplicates",
    Option  "v" ["verbose"]        (NoArg  set_verbose)       "Print more diagnostics",
    Option "h?" ["help","usage"]   (NoArg  (const usage))     "Print this message" ]

  where
    set_output   f c =                    return $ c { output = writeRawBamFile f } 
    set_qual     n c = readIO n >>= \a -> return $ c { collapse = cons_collapse' a }
    set_no_strand  c =                    return $ c { strand_preserved = False }
    set_verbose    c =                    return $ c { debug = hPutStr stderr }
    set_improper   c =                    return $ c { filter_enee = Just }
    set_single     c =                    return $ c { filter_enee = make_single }
    set_cheap      c =                    return $ c { collapse = cheap_collapse' }
    set_count_only c =                    return $ c { collapse = cheap_collapse', output = const skipToEof }
    set_keep       c =                    return $ c { keep_all = True }
    set_unaligned  c =                    return $ c { keep_unaligned = True }
    set_len      n c = readIO n >>= \a -> return $ c { min_len = a }
    set_no_rg      c =                    return $ c { get_label = get_no_library }

    set_range    a c
        | a == "A" || a == "a" = return $ c { which = All }
        | a == "U" || a == "u" = return $ c { which = Unaln }
        | otherwise = case reads a of
                [ (x,"")    ] -> return $ c { which = Some (Refseq $ x-1) (Refseq $ x-1) }
                [ (x,'-':b) ] -> readIO b >>= \y -> 
                                 return $ c { which = Some (Refseq $ x-1) (Refseq $ y-1) }
                _ -> fail $ "parse error in " ++ show a                                 

usage :: IO a
usage = do p <- getProgName
           hPutStrLn stderr $ "Usage: " ++ usageInfo (p ++ info) options 
           exitSuccess
  where 
    info = " [option...] [bam-file...]\n\
           \Removes PCR duplicates from BAM files and calls a consensus for each duplicate set.  \n\
           \Input files must be sorted by coordinate and are merged on the fly.  Options are:"
    
cons_collapse' :: Word8 -> Bool -> Collapse
cons_collapse' m False = cons_collapse m
cons_collapse' m True  = cons_collapse_keep m

cheap_collapse' :: Bool -> Collapse
cheap_collapse'  False = cheap_collapse
cheap_collapse'  True  = cheap_collapse_keep

-- | Get library from BAM record.
-- This gets the read group from a bam record, then the library for read
-- group.  This will work correctly if and only if the RG-LB field is
-- the name of the "Ur-Library", the common one before the first
-- amplification.
--
-- If no RG-LB field is present, RG-SM is used instead.  This will work
-- if and only if no libraries were aliquotted and then pooled again.
--
-- Else the RG-ID field is used.  This will work if and only if read
-- groups correspond directly to libraries.
--
-- If no RG is present, the empty string is returned.  This serves as
-- fall-back.

get_library, get_no_library :: M.Map Seqid Seqid -> BamRaw -> Seqid
get_library  tbl br = M.findWithDefault rg rg tbl where rg = br_extAsString "RG" br
get_no_library _  _ = S.empty

mk_rg_tbl :: BamMeta -> M.Map Seqid Seqid
mk_rg_tbl hdr = M.fromList
    [ (rg_id, rg_lb)
    | ('R','G',fields) <- meta_other_shit hdr
    , rg_id <- take 1   [ i | ('I','D',i) <- fields ]
    , rg_lb <- take 1 $ [ l | ('L','B',l) <- fields ]
                     ++ [ s | ('S','M',s) <- fields ]
                     ++ [ rg_id ] ]

data Counts = Counts { tin          :: !Int
                     , tout         :: !Int
                     , good_singles :: !Int
                     , good_total   :: !Int }

main :: IO ()
main = do
    args <- getArgs
    let (opts, files, errors) = getOpt Permute options args
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Conf{..} <- foldr (>=>) return opts defaults

    when (null args) $ do t <- hIsTerminalDevice stdout
                          when t $ hPutStrLn stderr "Cowardly refusing to write BAM to a terminal." >> exitFailure

    add_pg <- addPG $ Just version
    (counts, ()) <- mergeInputRanges which files >=> run $ \hdr -> do
       let tbl = mk_rg_tbl hdr
       unless (M.null tbl) $ liftIO $ do
                debug "mapping of read groups to libraries:\n"
                mapM_ debug [ unpackSeqid k ++ " --> " ++ unpackSeqid v ++ "\n" | (k,v) <- M.toList tbl ]

       output' <- takeWhileE is_halfway_aligned ><>
                  mapChunks (mapMaybe filter_enee . filter ((>= min_len) . eff_len)) ><>
                  progress debug (meta_refs hdr) ><>
                  rmdup (get_label tbl) strand_preserved (collapse keep_all) $ 
                  count_all (get_label tbl) `I.zip` output (add_pg hdr)
       if keep_unaligned
         then output'
         else lift (run output')

    hPutStr stderr . unlines $
        "\27[K#RG\tin\tout\tin@MQ20\tsingle@MQ20\tunseen\ttotal\t%unique\t%exhausted"
        : map (uncurry do_report) (M.toList counts)


do_report :: Seqid -> Counts -> String
do_report lbl Counts{..} = intercalate "\t" fs
  where
    fs = label : showNum tin : showNum tout : showNum good_total : showNum good_singles : 
         report_estimate (estimateComplexity good_total good_singles)

    label = if S.null lbl then "--" else unpackSeqid lbl

    report_estimate  Nothing                = [ "N/A" ]
    report_estimate (Just good_grand_total) =
            [ showOOM (grand_total - fromIntegral tout)
            , showOOM grand_total
            , showFFloat (Just 1) rate []
            , showFFloat (Just 1) exhaustion [] ]
      where 
        grand_total = good_grand_total * fromIntegral tout / fromIntegral good_total 
        exhaustion  = 100 * fromIntegral good_total / good_grand_total
        rate        = 100 * fromIntegral tout / fromIntegral tin :: Double


count_all :: (BamRaw -> Seqid) -> Iteratee [BamRaw] m (M.Map Seqid Counts)
count_all lbl = I.foldl' plus M.empty
  where
    plus m br = M.insert (lbl br) cs m
      where
        !cs = plus1 (M.findWithDefault (Counts 0 0 0 0) (lbl br) m) br

    plus1 (Counts ti to gs gt) br = Counts ti' to' gs' gt'
      where
        !ti' = ti + br_extAsInt 1 "XP" br
        !to' = to + 1
        !gs' = if br_mapq br >= 20 && br_extAsInt 1 "XP" br == 1 then gs + 1 else gs
        !gt' = if br_mapq br >= 20 then gt + 1 else gt


eff_len :: BamRaw -> Int
eff_len br | br_isProperlyPaired br = abs $ br_isize br
           | otherwise              = br_l_seq br

is_halfway_aligned :: BamRaw -> Bool
is_halfway_aligned br = not (br_isUnmapped br) || not (br_isMateUnmapped br)

is_aligned :: BamRaw -> Maybe BamRaw
is_aligned br | br_isUnmapped br || not (isValidRefseq (br_rname br)) = Nothing
              | br_isPaired br && br_isMateUnmapped br                = Nothing
              | otherwise                                             = Just br

make_single :: BamRaw -> Maybe BamRaw
make_single br | br_isPaired br && br_isSecondMate br = Nothing
               | br_isUnmapped br                     = Nothing
               | not (br_isPaired br)                 = Just br
               | otherwise = Just $! mutateBamRaw br $ do setFlag $ br_flag br .&. complement pair_flags
                                                          setMrnm $ invalidRefseq
                                                          setMpos $ invalidPos
                                                          setIsize $ 0
  where                                    
    pair_flags = flagPaired .|. flagProperlyPaired .|.
                 flagFirstMate .|. flagSecondMate .|.
                 flagMateUnmapped
                                

progress :: MonadIO m => (String -> IO ()) -> Refs -> Enumeratee [BamRaw] [BamRaw] m a
progress put refs = eneeCheckIfDone (liftI . go 0)
  where
    go !_ k (EOF         mx) = idone (liftI k) (EOF mx)
    go !n k (Chunk    [   ]) = liftI $ go n k
    go !n k (Chunk as@(a:_)) = do let !n' = n + length as
                                      nm = unpackSeqid (sq_name (getRef refs (br_rname a))) ++ ":"
                                  when (n `div` 65536 /= n' `div` 65536) $ liftIO $ put $
                                        "\27[KRmdup at " ++ nm ++ showNum (br_pos a) ++ "\r"
                                  eneeCheckIfDone (liftI . go n') . k $ Chunk as


mergeInputRanges :: MonadCatchIO m
    => Which -> [FilePath] -> Enumerator' BamMeta [BamRaw] m a
mergeInputRanges All      fps   = mergeInputs combineCoordinates fps
mergeInputRanges  _  [        ] = \k -> return $ k mempty
mergeInputRanges rng (fp0:fps0) = go fp0 fps0
  where
    enum1  fp k1 = case rng of All      -> decodeAnyBamFile       fp k1
                               Some x y -> decodeBamFileRange x y fp k1
                               -- XXX Unaln ->

    go fp [       ] = enum1 fp
    go fp (fp1:fps) = mergeEnums' (go fp1 fps) (enum1 fp) combineCoordinates

decodeBamFileRange :: MonadCatchIO m
                   => Refseq -> Refseq -> FilePath 
                   -> (BamMeta -> Iteratee [BamRaw] m a) 
                   -> m (Iteratee [BamRaw] m a)
decodeBamFileRange x y fp k0 = do
    idx <- liftIO $ readBamIndex fp
    enumFileRandom defaultBufSize fp >=> run $
        joinI $ decompressBgzf $ do
            hdr <- decodeBam return
            foldr ((>=>) . decodeBamSequence idx) return [x..y] $ hdr >>= k0

