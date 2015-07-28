{-# LANGUAGE RecordWildCards, BangPatterns, FlexibleContexts #-}
import Bio.Bam
import Bio.Base
import Bio.Util ( showNum, showOOM, estimateComplexity )
import Control.Monad
import Control.Monad.ST ( runST )
import Data.Bits
import Data.Foldable ( toList )
import Data.List ( intercalate )
import Data.Maybe
import Data.Monoid ( mempty )
import Data.Ord ( comparing )
import Data.Vector.Algorithms.Intro ( sortBy )
import Data.Version ( showVersion )
import Numeric ( showFFloat )
import Paths_biohazard ( version )
import System.Console.GetOpt
import System.Environment ( getArgs, getProgName )
import System.Exit
import System.IO

import qualified Data.ByteString.Char8  as S
import qualified Data.HashMap.Strict    as M
import qualified Data.IntMap            as IM
import qualified Data.Iteratee          as I
import qualified Data.Sequence          as Z
import qualified Data.Vector            as V

data Conf = Conf {
    output :: Maybe ((BamRaw -> Seqid) -> BamMeta -> Iteratee [BamRaw] IO ()),
    strand_preserved :: Bool,
    collapse :: Bool -> Collapse,
    clean_multimap :: BamRaw -> IO (Maybe BamRaw),
    keep_all :: Bool,
    keep_unaligned :: Bool,
    keep_improper :: Bool,
    transform :: BamRaw -> Maybe BamRaw,
    min_len :: Int,
    min_qual :: Qual,
    get_label :: M.HashMap Seqid Seqid -> BamRaw -> Seqid,
    putResult :: String -> IO (),
    debug :: String -> IO (),
    which :: Which,
    circulars :: Refs -> IO (IM.IntMap (Seqid,Int), Refs) }

-- | Which reference sequences to scan
data Which = All | Some Refseq Refseq | Unaln deriving Show

defaults :: Conf
defaults = Conf { output = Nothing
                , strand_preserved = True
                , collapse = cons_collapse' (Q 60)
                , clean_multimap = check_flags
                , keep_all = False
                , keep_unaligned = False
                , keep_improper = False
                , transform = Just
                , min_len = 0
                , min_qual = Q 0
                , get_label = get_library
                , putResult = putStr
                , debug = \_ -> return ()
                , which = All
                , circulars = \rs -> return (IM.empty, rs) }

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option  "o" ["output"]         (ReqArg set_output "FILE") "Write to FILE (default: no output, count only)",
    Option  "O" ["output-lib"]     (ReqArg set_lib_out "PAT") "Write each lib to file named following PAT",
    Option  [ ] ["debug"]          (NoArg  set_debug_out)     "Write textual debugging output",
    Option  "z" ["circular"]       (ReqArg add_circular "CHR:LEN") "Refseq CHR is circular with length LEN",
    Option  "R" ["refseq"]         (ReqArg set_range "RANGE") "Read only range of reference sequences",
    Option  "p" ["improper-pairs"] (NoArg  set_improper)      "Include improper pairs",
    Option  "u" ["unaligned"]      (NoArg  set_unaligned)     "Include unaligned reads and pairs",
    Option  "1" ["single-read"]    (NoArg  set_single)        "Pretend there is no second mate",
    Option  "m" ["multimappers"]   (NoArg  set_multi)         "Process multi-mappers (by dropping secondary alignments)",
    Option  "c" ["cheap"]          (NoArg  set_cheap)         "Cheap computation: skip the consensus calling",
    Option  "k" ["keep","mark-only"](NoArg set_keep)          "Mark duplicates, but include them in output",
    Option  "Q" ["max-qual"]       (ReqArg set_qual "QUAL")   "Set maximum quality after consensus call to QUAL",
    Option  "l" ["min-length"]     (ReqArg set_len "LEN")     "Discard reads shorter than LEN",
    Option  "q" ["min-mapq"]       (ReqArg set_mapq "QUAL")   "Discard reads with map quality lower than QUAL",
    Option  "s" ["no-strand"]      (NoArg  set_no_strand)     "Strand of alignments is uninformative",
    Option  "r" ["ignore-rg"]      (NoArg  set_no_rg)         "Ignore read groups when looking for duplicates",
    Option  "v" ["verbose"]        (NoArg  set_verbose)       "Print more diagnostics",
    Option "h?" ["help","usage"]   (NoArg  (const usage))     "Display this message",
    Option "V"  ["version"]        (NoArg  (const vrsn))      "Display version number and exit" ]

  where
    set_output "-" c =                    return $ c { output = Just $ \_ -> pipeRawBamOutput, putResult = hPutStr stderr }
    set_output   f c =                    return $ c { output = Just $ \_ -> writeRawBamFile f }
    set_lib_out  f c =                    return $ c { output = Just $       writeLibBamFiles f }
    set_debug_out  c =                    return $ c { output = Just $ \_ -> pipeRawSamOutput, putResult = hPutStr stderr }
    set_qual     n c = readIO n >>= \q -> return $ c { collapse = cons_collapse' (Q q) }
    set_no_strand  c =                    return $ c { strand_preserved = False }
    set_verbose    c =                    return $ c { debug = hPutStr stderr }
    set_improper   c =                    return $ c { keep_improper = True }
    set_single     c =                    return $ c { transform = make_single }
    set_cheap      c =                    return $ c { collapse = cheap_collapse' }
    set_keep       c =                    return $ c { keep_all = True }
    set_unaligned  c =                    return $ c { keep_unaligned = True }
    set_len      n c = readIO n >>= \l -> return $ c { min_len = l }
    set_mapq     n c = readIO n >>= \q -> return $ c { min_qual = Q q }
    set_no_rg      c =                    return $ c { get_label = get_no_library }
    set_multi      c =                    return $ c { clean_multimap = clean_multi_flags }

    set_range    a c
        | a == "A" || a == "a" = return $ c { which = All }
        | a == "U" || a == "u" = return $ c { which = Unaln }
        | otherwise = case reads a of
                [ (x,"")    ] -> return $ c { which = Some (Refseq $ x-1) (Refseq $ x-1) }
                [ (x,'-':b) ] -> readIO b >>= \y ->
                                 return $ c { which = Some (Refseq $ x-1) (Refseq $ y-1) }
                _ -> fail $ "parse error in " ++ show a

    add_circular a c = case break ((==) ':') a of
        (nm,':':r) -> case reads r of
            [(l,[])] | l > 0 -> return $ c { circulars = add_circular' (S.pack nm) l (circulars c) }
            _ -> fail $ "couldn't parse length " ++ show r ++ " for " ++ show nm
        _ -> fail $ "couldn't parse \"circular\" argument " ++ show a

    add_circular' nm l io refs = do
        (m1, refs') <- io refs
        case filter (S.isPrefixOf nm . sq_name . snd) $ zip [0..] $ toList refs' of
            [(k,a)] | sq_length a >= l -> let m2     = IM.insert k (sq_name a,l) m1
                                              refs'' = Z.update k (a { sq_length = l }) refs'
                                          in return (m2, refs'')
                    | otherwise -> fail $ "cannot wrap " ++ show nm ++ " to " ++ show l
                                       ++ ", which is more than the original " ++ show (sq_length a)
            [] -> fail $ "no match for target sequence " ++ show nm
            _ -> fail $ "target sequence " ++ show nm ++ " is ambiguous"

vrsn :: IO a
vrsn = do pn <- getProgName
          hPutStrLn stderr $ pn ++ ", version " ++ showVersion version
          exitSuccess

usage :: IO a
usage = do p <- getProgName
           hPutStrLn stderr $ "Usage: " ++ usageInfo (p ++ info) options
           exitSuccess
  where
    info = " [option...] [bam-file...]\n\
           \Removes PCR duplicates from BAM files and calls a consensus for each duplicate set.  \n\
           \Input files must be sorted by coordinate and are merged on the fly.  Options are:"

cons_collapse' :: Qual -> Bool -> Collapse
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

get_library, get_no_library :: M.HashMap Seqid Seqid -> BamRaw -> Seqid
get_library  tbl br = M.lookupDefault rg rg tbl where rg = br_extAsString "RG" br
get_no_library _  _ = S.empty

mk_rg_tbl :: BamMeta -> M.HashMap Seqid Seqid
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
    when (null args) usage
    let (opts, files, errors) = getOpt Permute options args
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Conf{..} <- foldr (>=>) return opts defaults

    add_pg <- addPG $ Just version
    (counts, ()) <- mergeInputRanges which files >=> run $ \hdr -> do
       (circtable, refs') <- liftIO $ circulars (meta_refs hdr)
       let tbl = mk_rg_tbl hdr
       unless (M.null tbl) $ liftIO $ do
                debug "mapping of read groups to libraries:\n"
                mapM_ debug [ unpackSeqid k ++ " --> " ++ unpackSeqid v ++ "\n" | (k,v) <- M.toList tbl ]

       let filters = mapChunks (mapMaybe transform) ><>
                     mapChunksM (mapMM clean_multimap) ><>
                     filterStream (\br -> (keep_unaligned || is_aligned br) &&
                                          (keep_improper || is_proper br) &&
                                          eff_len br >= min_len) ><>
                     progressPos "Rmdup at " debug refs'

       let (co, ou) = case output of Nothing -> (cheap_collapse', skipToEof)
                                     Just  o -> (collapse, joinI $ wrapSortWith circtable $
                                                           o (get_label tbl) (add_pg hdr { meta_refs = refs' }))

       ou' <- takeWhileE is_halfway_aligned ><> filters ><>
              normalizeSortWith circtable ><>
              filterStream (\br -> br_mapq br >= min_qual) ><>
              rmdup (get_label tbl) strand_preserved (co keep_all) $
              count_all (get_label tbl) `I.zip` ou

       liftIO $ debug "\27[Krmdup done; copying junk\n"

       case which of
            Unaln              -> joinI $ filters ou'
            _ | keep_unaligned -> joinI $ filters ou'
            _                  -> lift (run ou')

    putResult . unlines $
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


-- | Counting reads:  we count total read in (ti), total reads out (to),
-- good (confidently mapped) singletons out (gs), total good
-- (confidently mapped) reads out (gt).  Paired reads count 1, unpaired
-- reads count 2, and at the end we divide by 2.  This ensures that we
-- don't double count mate pairs, while still working mostly sensibly in
-- the presence of broken BAM files.

count_all :: Functor m => (BamRaw -> Seqid) -> Iteratee [BamRaw] m (M.HashMap Seqid Counts)
count_all lbl = M.map fixup `fmap` I.foldl' plus M.empty
  where
    plus m br = M.insert (lbl br) cs m
      where
        !cs = plus1 (M.lookupDefault (Counts 0 0 0 0) (lbl br) m) br

    plus1 (Counts ti to gs gt) br = Counts ti' to' gs' gt'
      where
        !w   = if br_isPaired br then 1 else 2
        !ti' = ti + w * br_extAsInt 1 "XP" br
        !to' = to + w
        !gs' = if br_mapq br >= Q 20 && br_extAsInt 1 "XP" br == 1 then gs + w else gs
        !gt' = if br_mapq br >= Q 20                               then gt + w else gt

    fixup (Counts ti to gs gt) = Counts (div ti 2) (div to 2) (div gs 2) (div gt 2)

eff_len :: BamRaw -> Int
eff_len br | br_isProperlyPaired br = abs $ br_isize br
           | otherwise              = br_l_seq br

is_halfway_aligned :: BamRaw -> Bool
is_halfway_aligned br = not (br_isUnmapped br) || not (br_isMateUnmapped br)

is_aligned :: BamRaw -> Bool
is_aligned br = not (br_isUnmapped br && br_isMateUnmapped br) && isValidRefseq (br_rname br)

is_proper :: BamRaw -> Bool
is_proper br = not (br_isPaired br) ||
               (br_isMateUnmapped br == br_isUnmapped br && br_isProperlyPaired br)

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


mergeInputRanges :: (MonadIO m, MonadMask m)
                 => Which -> [FilePath] -> Enumerator' BamMeta [BamRaw] m a
mergeInputRanges All      fps   = mergeInputs combineCoordinates fps
mergeInputRanges  _  [        ] = \k -> return $ k mempty
mergeInputRanges rng (fp0:fps0) = go fp0 fps0
  where
    enum1  fp k1 = case rng of All      -> decodeAnyBamFile                 fp k1
                               Some x y -> decodeBamFileRange           x y fp k1
                               Unaln    -> decodeWithIndex eneeBamUnaligned fp k1

    go fp [       ] = enum1 fp
    go fp (fp1:fps) = mergeEnums' (go fp1 fps) (enum1 fp) combineCoordinates

    decodeBamFileRange x y = decodeWithIndex $
            \idx -> foldr ((>=>) . eneeBamRefseq idx) return [x..y]


decodeWithIndex :: (MonadIO m, MonadMask m)
                => (BamIndex () -> Enumeratee [BamRaw] [BamRaw] m a)
                -> FilePath -> (BamMeta -> Iteratee [BamRaw] m a)
                -> m (Iteratee [BamRaw] m a)
decodeWithIndex enum fp k0 = do
    idx <- liftIO $ readBamIndex fp
    decodeAnyBamFile fp >=> run $ enum idx . k0


writeLibBamFiles :: (MonadIO m, MonadMask m)
                 => FilePath -> (BamRaw -> Seqid) -> BamMeta -> Iteratee [BamRaw] m ()
writeLibBamFiles fp lbl hdr = tryHead >>= loop M.empty
  where
    loop m  Nothing  = liftIO . mapM_ run $ M.elems m
    loop m (Just br) = do
        let !l = lbl br
        let !it = M.lookupDefault (writeRawBamFile (fp `subst` l) hdr) l m
        it' <- liftIO $ enumPure1Chunk [br] it
        let !m' = M.insert l it' m
        tryHead >>= loop m'

    subst [            ] _ = []
    subst ('%':'s':rest) l = unpackSeqid l ++ subst rest l
    subst ('%':'%':rest) l = '%' : subst rest l
    subst ( c :    rest) l =  c  : subst rest l


mapMM :: Monad m => (a -> m (Maybe b)) -> [a] -> m [b]
mapMM f = go []
  where
    go acc [    ] = return $ reverse acc
    go acc (a:as) = do b <- f a ; go (maybe acc (:acc) b) as


check_flags :: Monad m => BamRaw -> m (Maybe BamRaw)
check_flags b | br_extAsInt 1 "HI" b /= 1 = fail "cannot deal with HI /= 1"
              | br_extAsInt 1 "IH" b /= 1 = fail "cannot deal with IH /= 1"
              | br_extAsInt 1 "NH" b /= 1 = fail "cannot deal with NH /= 1"
              | otherwise                 = return $ Just b

clean_multi_flags :: Monad m => BamRaw -> m (Maybe BamRaw)
clean_multi_flags b = return $ if br_extAsInt 1 "HI" b /= 1 then Nothing else Just b'
  where
    b' = mutateBamRaw b $ mapM_ removeExt ["HI","IH","NH"]


-- Given a map from reference sequences to arguments, extract those
-- groups as list, apply a function to the argument and the list, pass
-- the result on.  Absent groups are passed on as they are.  Note that
-- ordering within groups is messed up (it doesn't matter here).
mapAtGroups :: Monad m => IM.IntMap a -> (a -> [BamRaw] -> [BamRaw]) -> Enumeratee [BamRaw] [BamRaw] m b
mapAtGroups m f = eneeCheckIfDonePass no_group
  where
    no_group k (Just e) = idone (liftI k) $ EOF (Just e)
    no_group k Nothing  = tryHead >>= maybe (idone (liftI k) $ EOF Nothing) (\a -> no_group_1 a k Nothing)

    no_group_1 _ k (Just e) = idone (liftI k) $ EOF (Just e)
    no_group_1 a k Nothing  = case IM.lookup (br_rname_int a) m of
            Nothing  -> eneeCheckIfDonePass no_group . k $ Chunk [a]
            Just arg -> cont_group (br_rname a) arg [a] k Nothing

    cont_group rn arg acc k (Just e) = idone (liftI k) $ EOF (Just e)
    cont_group rn arg acc k Nothing = tryHead >>= maybe flush_eof check1
      where
        flush_eof  = idone (k $ Chunk $ f arg acc) (EOF Nothing)
        flush_go a = eneeCheckIfDonePass (no_group_1 a) . k . Chunk $ f arg acc
        check1 a | br_rname a == rn = cont_group rn arg (a:acc) k Nothing
                 | otherwise        = flush_go a

    br_rname_int = fromIntegral . unRefseq . br_rname

normalizeSortWith :: Monad m => IM.IntMap (Seqid, Int) -> Enumeratee [BamRaw] [BamRaw] m a
normalizeSortWith m = mapAtGroups m $ \(nm,l) -> sortPos . map (normalizeTo nm l)

wrapSortWith :: Monad m => IM.IntMap (Seqid, Int) -> Enumeratee [BamRaw] [BamRaw] m a
wrapSortWith m = mapAtGroups m $ \(_,l) -> sortPos . concatMap (wrapTo l)

sortPos :: [BamRaw] -> [BamRaw]
sortPos l = V.toList $ runST (V.unsafeThaw (V.fromList l) >>= \vm -> sortBy (comparing br_pos) vm >> V.unsafeFreeze vm)
