{-# LANGUAGE RecordWildCards, BangPatterns #-}
import Bio.File.Bam
import Bio.File.Bam.Rmdup
import Bio.File.Bam.Fastq ( removeWarts )
import Bio.Iteratee
import Control.Monad
import System.Console.GetOpt
import System.Environment ( getArgs, getProgName )
import System.Exit
import System.IO

import qualified Data.Map as M

data Conf = Conf {
    output :: BamMeta -> Iteratee [BamRec] IO (),
    max_qual :: Int,
    strand_preserved :: Bool,
    debug :: String -> IO () }

defaults :: Conf
defaults = Conf { output = writeBamHandle stdout
                , max_qual = 60
                , strand_preserved = True
                , debug = \_ -> return () }

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option  "o" ["output"]       (ReqArg set_output "FILE") "Write to FILE (default: stdout)",
    Option  "Q" ["max-qual"]     (ReqArg set_qual "QUAL")   "Set maximum quality after consensus call to QUAL",
    Option  "s" ["no-strand"]    (NoArg  set_no_strand)     "Strand of alignments is uninformative",
    Option  "v" ["verbose"]      (NoArg  set_verbose)       "Print more diagnostics",
    Option "h?" ["help","usage"] (NoArg  usage)             "Print this message" ]
  where
    set_output  f c =                    return $ c { output = writeBamFile f } 
    set_qual    n c = readIO n >>= \a -> return $ c { max_qual = a }
    set_no_strand c =                    return $ c { strand_preserved = False }
    set_verbose   c =                    return $ c { debug = hPutStr stderr }

    usage _ = do p <- getProgName
                 hPutStrLn stderr $ usageInfo (p ++ info)  options 
                 exitSuccess
    info = " [option...] [bam-file...]\n\
           \Removes PCR duplicates from BAM files and calls a consensus for each duplicate set.  \
           \Input files must be sorted by coordinate and are merged on the fly.  \
           \Options are:"
    

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

get_library :: M.Map Seqid Seqid -> BamRec -> Seqid
get_library tbl = \br -> let rg = extAsString "RG" br in M.findWithDefault rg rg tbl

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
    (opts, files, errors) <- getOpt Permute options `fmap` getArgs
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Conf{..} <- foldr (>=>) return opts defaults
    enum_all_input_files files >=> run $ \hdr -> do
       let tbl = mk_rg_tbl hdr
       unless (M.null tbl) $ liftIO $ do
                debug "mapping of read groups to libraries:\n"
                mapM_ debug [ unpackSeqid k ++ " --> " ++ unpackSeqid v ++ "\n" | (k,v) <- M.toList tbl ]

       joinI $ takeWhileE is_halfway_aligned $
           joinI $ filterStream is_aligned $ 
           joinI $ progress debug (meta_refs hdr) $
           joinI $ rmdup (get_library tbl) strand_preserved (fromIntegral $ min 255 max_qual) $
           output hdr

is_halfway_aligned :: BamRec -> Bool
is_halfway_aligned br = not (isUnmapped br) || not (isMateUnmapped br)

is_aligned :: BamRec -> Bool
is_aligned br | not (isValidRefseq (b_rname br)) = False 
              | isPaired br = not (isUnmapped br) && not (isMateUnmapped br)
              | otherwise   = not (isUnmapped br)


enum_all_input_files :: [FilePath] -> Enumerator' BamMeta [BamRec] IO a
enum_all_input_files [        ] = enum_input_file "-"
enum_all_input_files (fp0:fps0) = go fp0 fps0
  where
    go fp [       ] = enum_input_file fp
    go fp (fp1:fps) = go fp1 fps ? enum_input_file fp 
    a ? b = mergeEnums' a b (const combineCoordinates)

basicFilters :: Monad m => Enumeratee [BamRec] [BamRec] m a
basicFilters = mapStream removeWarts ><> filterStream checkFlags

checkFlags :: BamRec -> Bool
checkFlags br | isPaired br = not $ isUnmapped br && isMateUnmapped br
              | otherwise   = not $ isFailsQC br || isUnmapped br


enum_input_file :: MonadCatchIO m => FilePath -> Enumerator' BamMeta [BamRec] m a
enum_input_file f it = enum_input_file' f >=> run $ \hdr -> basicFilters $ it hdr

enum_input_file' :: MonadCatchIO m => FilePath -> Enumerator' BamMeta [BamRec] m a
enum_input_file' "-"  = (enumHandle defaultBufSize stdin >=> run) . decodeAnyBamOrSam
enum_input_file' path = decodeAnyBamOrSamFile path


progress :: MonadIO m => (String -> IO ()) -> Refs -> Enumeratee [BamRec] [BamRec] m a
progress put refs = eneeCheckIfDone (liftI . go 0)
  where
    go !_ k (EOF mx) = do liftIO $ put "\27[KDone.\n" 
                          idone (liftI k) (EOF mx)
    
    go !n k (Chunk []) = liftI $ go n k
    go !n k (Chunk as@(a:_)) = do let !n' = n + length as
                                      nm = shows (sq_name (getRef refs (b_rname a))) ":"
                                  when (n `div` 16384 /= n' `div` 16384) $ liftIO $ put $
                                        "\27[KRmdup at " ++ nm ++ shows (b_pos a) "\r"
                                  eneeCheckIfDone (liftI . go n') . k $ Chunk as

