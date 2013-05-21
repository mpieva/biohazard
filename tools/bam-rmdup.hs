{-# LANGUAGE RecordWildCards, BangPatterns #-}
import Bio.File.Bam
import Bio.File.Bam.Rmdup
import Bio.Iteratee
import Bio.Util ( showNum, showOOM, estimateComplexity )
import Control.Monad
import Data.Bits
import Data.Maybe
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
    filter_enee :: BamRaw -> Maybe BamRaw,
    min_len :: Int,
    get_label :: M.Map Seqid Seqid -> BamRaw -> Seqid,
    debug :: String -> IO () }

defaults :: Conf
defaults = Conf { output = pipeRawBamOutput
                , strand_preserved = True
                , collapse = cons_collapse' 60
                , keep_all = False
                , filter_enee = is_aligned
                , min_len = 0
                , get_label = get_library
                , debug = \_ -> return () }

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option  "o" ["output"]         (ReqArg set_output "FILE") "Write to FILE (default: stdout)",
    Option  "p" ["improper-pairs"] (NoArg  set_improper)      "Include improper pairs",
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
    set_len      n c = readIO n >>= \a -> return $ c { min_len = a }
    set_no_rg      c =                    return $ c { get_label = get_no_library }


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
get_library  tbl = \br -> let rg = br_extAsString "RG" br in M.findWithDefault rg rg tbl
get_no_library _ = \_  -> S.empty

mk_rg_tbl :: BamMeta -> M.Map Seqid Seqid
mk_rg_tbl hdr = M.fromList
    [ (rg_id, rg_lb)
    | ('R','G',fields) <- meta_other_shit hdr
    , rg_id <- take 1 [ i | ('I','D',i) <- fields ]
    , rg_lb <- take 1 $ [ l | ('L','B',l) <- fields ]
                     ++ [ s | ('S','M',s) <- fields ]
                     ++ [ rg_id ] ]

main :: IO ()
main = do
    args <- getArgs
    let (opts, files, errors) = getOpt Permute options args
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Conf{..} <- foldr (>=>) return opts defaults

    when (null args) $ do t <- hIsTerminalDevice stdout
                          when t $ hPutStrLn stderr "Cowardly refusing to write BAM to a terminal." >> usage

    add_pg <- addPG $ Just version
    (tin, (tout, good_singles, good_total, ())) <- mergeInputs combineCoordinates files >=> run $ \hdr -> do
       let tbl = mk_rg_tbl hdr
       unless (M.null tbl) $ liftIO $ do
                debug "mapping of read groups to libraries:\n"
                mapM_ debug [ unpackSeqid k ++ " --> " ++ unpackSeqid v ++ "\n" | (k,v) <- M.toList tbl ]

       joinI $ takeWhileE is_halfway_aligned $
           joinI $ mapChunks (mapMaybe filter_enee . filter ((>= min_len) . eff_len)) $
           joinI $ progress debug (meta_refs hdr) $
           I.zip I.length $ joinI $ rmdup (get_label tbl) strand_preserved (collapse keep_all) $
                            I.zip4 I.length count_singles count_good $ output (add_pg hdr)


        
    let report_estimate Nothing = "Complexity cannot be estimated.\n"
        report_estimate (Just good_grand_total) =
            "Estimate " ++ showOOM (grand_total - fromIntegral tout) ++ " dark matter molecules, " ++
            showOOM grand_total ++ " total distinct molecules.\n" ++
            "Unique fraction so far is " ++ showFFloat (Just 1) rate 
            "%, library is " ++ showFFloat (Just 1) exhaustion "% exhausted.\n"
          where 
            grand_total = good_grand_total * fromIntegral tout / fromIntegral good_total 
            exhaustion  = 100 * fromIntegral good_total / good_grand_total
            rate        = 100 * fromIntegral tout / fromIntegral tin                        :: Double

    hPutStr stderr $ "\27[KDone; " ++ showNum (tin::Int) ++
                     " aligned reads in, " ++ showNum (tout::Int) ++ 
                     " aligned reads out; " ++ showNum good_total ++
                     " at MAPQ>=20, " ++ showNum good_singles ++
                     " singletons at MAPQ>=20.\n" ++
                     report_estimate (estimateComplexity good_total good_singles)


count_if :: Monad m => (a -> Bool) -> Iteratee [a] m Int
count_if p = I.foldl' (\acc a -> if p a then acc + 1 else acc) 0

count_singles :: Monad m => Iteratee [BamRaw] m Int
count_singles = count_if $ \br -> br_mapq br >= 20 && br_extAsInt 1 "XP" br == 1

count_good :: Monad m => Iteratee [BamRaw] m Int
count_good = count_if $ \br -> br_mapq br >= 20


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
