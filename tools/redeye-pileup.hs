{-# LANGUAGE RecordWildCards, BangPatterns, OverloadedStrings, FlexibleContexts #-}
-- Command line driver for simple genotype calling.  We have three
-- separate steps:  Pileup from a BAM file (or multiple merged files) to
-- produce likelihoods (and some auxillary statistics).  These are
-- written into an Avro container.  Next we need to estimate parameters,
-- in the simplest case divergence and heterozygosity.  We can save some
-- time by fusing this with the first step.  The final step is calling
-- bases by scnaning the Avro container and applying some model, and
-- again, in the simplest case that's just divergence and
-- heterozygosity.  We keep that separate, because different models will
-- require different programs.  So here we produce likelihoods and
-- a simple model fit.

-- The likelihoods depend on damage parameters and an error model,
-- otherwise they are 'eternal'.  (For the time being, it's probably
-- wise to go with the naïve error model.)  Technically, they also
-- depend on ploidy, but since only diploid organisms are interesting
-- right now, we fix that to two.  We pay some overhead on the sex
-- chromosomes, but the simplification is worth it.

-- About damage parameters:  We effectively have three different models
-- (SS, DS, no damage) and it may not be possible to choose one a
-- priori.  To manage this cleanly, we should have one universal model,
-- but the three we have are not generalizations of each other.
-- However, all can be generalized into one model with slightly more
-- parameters.  See tools/dmg-est.hs for how we fit the model.

-- Calling is always diploid, for maximum flexibility.  We don't really
-- support higher ploidies, so the worst damage is that we output an
-- overhead of 150% useless likelihood values for the sex chromosomes
-- and maybe estimate heterozygosity where there is none.

import Bio.Adna
import Bio.Bam
import Bio.Bam.Pileup
import Bio.Genocall
import Bio.Genocall.AvroFile
import Bio.Genocall.Metadata
import Bio.Prelude
import Data.Aeson
import Data.Avro
import Data.Vec.Packed               ( packMat )
import System.Console.GetOpt
import System.Directory              ( renameFile )
import System.FilePath
import System.IO

import qualified Data.ByteString.Char8          as S
import qualified Data.ByteString.Lazy           as BL
import qualified Data.Foldable                  as F
import qualified Data.HashMap.Strict            as H
import qualified Data.Text                      as T
import qualified Data.Text.Encoding             as T
import qualified Data.Vector.Storable           as VS
import qualified Data.Vector                    as V
import qualified Data.Vector.Unboxed            as U
import qualified Data.Vector.Unboxed.Mutable    as M
import qualified Data.Sequence                  as Z

data Conf = Conf {
    -- Generator for output file name.  Receives sample name and
    -- (optional) region as arguments.
    conf_output      :: String -> Maybe String -> FilePath,
    conf_metadata    :: FilePath,
    conf_theta       :: Maybe Double,
    conf_report      :: String -> IO () }

defaultConf :: Conf
defaultConf = Conf default_out (error "no metadata file specified") Nothing (\_ -> return ())
  where
    default_out smp   Nothing  = smp <> ".av"
    default_out smp (Just rgn) = smp <> "-" <> rgn <> ".av"

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "c"  ["config"] (ReqArg set_conf   "FILE") "Set name of json config file to FILE",
    Option "o"  ["output"] (ReqArg set_output "FILE") "Set out file schema to FILE",
    Option "t"  dep_param  (ReqArg set_theta  "FRAC") "Set dependency coefficient to FRAC (\"N\" to turn off)",
    Option "v"  ["verbose"]        (NoArg be_verbose) "Print more diagnostics",
    Option "h?" ["help","usage"]   (NoArg disp_usage) "Display this message" ]
  where
    dep_param = ["theta","dependency-coefficient"]

    disp_usage    _ = do pn <- getProgName
                         let blah = "Usage: " ++ pn ++ " [OPTION...] [SAMPLE [REGION...] ...]"
                         putStrLn $ usageInfo blah options
                         exitSuccess

    be_verbose    c = return $ c { conf_report = hPutStrLn stderr }
    set_conf   fn c = return $ c { conf_metadata = fn }

    set_theta "N" c = return $ c { conf_theta  = Nothing }
    set_theta   a c = (\t -> c { conf_theta       = Just   t }) <$> readIO a

    set_output fn c = return $ c { conf_output = mkoutput fn }

mkoutput :: FilePath -> String -> Maybe String -> FilePath
mkoutput str smp rgn = go str
  where
    go ('%':'s':s) = smp ++ go s
    go ('%':'r':s) = maybe id (++) rgn $ go s
    go ('%':'%':s) = '%' : go s
    go ('%': c :s) =  c  : go s
    go (     c :s) =  c  : go s
    go [         ] = [ ]

main :: IO ()
main = do
    (opts, samprgns, errs) <- getOpt Permute options <$> getArgs
    Conf{..} <- foldl (>>=) (return defaultConf) opts
    unless (null errs) $ mapM_ (hPutStrLn stderr) errs >> exitFailure

    -- samprgns contains samples and regions.  We define anything found in the metadata as sample,
    -- anything else as region to apply to the previous sample.  Name a sample "chr1" and you get
    -- what you deserve.

    samples <- flip split_sam_rgns samprgns <$> readMetadata conf_metadata
    when (null samples) $ hPutStrLn stderr "need (at least) one sample name" >> exitFailure

    forM_ samples $ \(sample, rgns) -> do
        meta <- readMetadata conf_metadata

        case H.lookup (fromString sample) meta of
            Nothing -> hPutStrLn stderr $ "unknown sample " ++ show sample

            Just smp -> forM_ rgns $ \rgn -> do
                let outstem = conf_output sample rgn
                    outfile = takeDirectory conf_metadata </> outstem
                    tmpfile = outfile ++ ".#"
                (tab,()) <- withFile tmpfile WriteMode                                                      $ \ohdl ->
                            mergeLibraries conf_report conf_metadata (sample_libraries smp) rgn >=> run     $ \hdr ->
                            progressPos (\(rs, p, _) -> (rs, p)) "GT call at " conf_report (meta_refs hdr) =$
                            pileup                                                                         =$
                            mapStream (calls conf_theta)                                                   =$
                            zipStreams tabulateSingle (output_avro ohdl $ meta_refs hdr)

                let upd_sample s = s { sample_div_tables = H.insert (maybe T.empty T.pack rgn)                 tab  (sample_div_tables s)
                                     , sample_avro_files = H.insert (maybe T.empty T.pack rgn) (fromString outstem) (sample_avro_files s) }

                updateMetadata (H.adjust upd_sample (fromString sample)) conf_metadata
                renameFile tmpfile outfile

mergeLibraries :: (MonadIO m, MonadMask m)
               => (String -> IO ()) -> FilePath
               -> [Library] -> Maybe String -> Enumerator' BamMeta [PosPrimChunks] m b
mergeLibraries       _   _ [    ]    _ = error "Sample must have at least one library."
mergeLibraries  report cfg [ l  ] mrgn = enumLibrary report cfg l mrgn
mergeLibraries  report cfg (l:ls) mrgn = mergeEnums' (mergeLibraries report cfg ls mrgn) (enumLibrary report cfg l mrgn) mm
  where
    mm _ = mergeSortStreams $ \(rs1, p1, _) (rs2, p2, _) -> if (rs1, p1) < (rs2, p2) then Less else NotLess

enumLibrary :: (MonadIO m, MonadMask m)
            => (String -> IO ()) -> FilePath
            -> Library -> Maybe String -> Enumerator' BamMeta [PosPrimChunks] m b
enumLibrary report cfg (Library nm fs mdp) mrgn output = do
    let (msg, dmg) = case mdp of Nothing -> ("no damage model", noDamage)
                                 Just dp -> ("universal damage parameters" ++ show dp, univDamage dp)

    liftIO . report $ "using " ++ msg ++ " for " ++ unpack nm
    mergeInputRgns mrgn combineCoordinates (map ((</>) (takeDirectory cfg) . unpack) fs)
        $== takeWhileE (isValidRefseq . b_rname . unpackBam)
        $== mapMaybeStream (\br ->
                let b = unpackBam br
                    m = dmg (isReversed b) (VS.length (b_qual b))
                in decompose (map packMat $ V.toList m) br)
        $ output

mergeInputRgns :: (MonadIO m, MonadMask m)
    => Maybe String
    -> (BamMeta -> Enumeratee [BamRaw] [BamRaw] (Iteratee [BamRaw] m) a)
    -> [FilePath] -> Enumerator' BamMeta [BamRaw] m a
mergeInputRgns     _      _  [        ] = \k -> return (k mempty)
mergeInputRgns  Nothing  (?)      fps   = mergeInputs (?) fps
mergeInputRgns (Just rs) (?) (fp0:fps0) = go fp0 fps0
  where
    enum1  fp k1 = do idx <- liftIO $ readBamIndex fp
                      enumFileRandom defaultBufSize fp >=> run >=> run $
                            decodeAnyBam $ \hdr ->
                            let Just ri = Z.findIndexL ((==) rs . unpack . sq_name) (meta_refs hdr)
                            in eneeBamRefseq idx (Refseq $ fromIntegral ri) $ k1 hdr

    go fp [       ] = enum1 fp
    go fp (fp1:fps) = mergeEnums' (go fp1 fps) (enum1 fp) (?)


-- | Ploidy is hardcoded as two here.  Can be changed if the need
-- arises.
--
-- XXX  For the time being, forward and reverse piles get concatenated.
-- For the naive call, this doesn't matter.  For the MAQ call, it feels
-- more correct to treat them separately and multiply (add?) the results.

calls :: Maybe Double -> Pile -> Calls
calls Nothing pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !s = simple_snp_call fq 2 $ uncurry (++) $ p_snp_pile pile
    !i = simple_indel_call 2 $ p_indel_pile pile
    -- XXX this should be a cmdline option
    -- fq = min 1 . (*) 1.333 . fromQual
    fq = fromQual

calls (Just theta) pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !i = simple_indel_call 2 $ p_indel_pile pile

    -- This lumps the two strands together
    -- !s = maq_snp_call 2 theta $ uncurry (++) $ p_snp_pile pile -- XXX

    -- This treats them separately
    !s | r == r'    = Snp_GLs (U.zipWith (*) x y) r     -- same ref base (normal case): multiply
       | r == nucsN = Snp_GLs y r'                      -- forward ref is N, use backward call
       | otherwise  = Snp_GLs x r                       -- else use forward call (even if this is incorrect,
      where                                             -- there is nothing else we can do here)
        Snp_GLs x r  = maq_snp_call 2 theta $ fst $ p_snp_pile pile
        Snp_GLs y r' = maq_snp_call 2 theta $ snd $ p_snp_pile pile


-- | Serialize the results from genotype calling in a sensible way.  We
-- write an Avro file, but we add another blocking layer on top so we
-- don't need to endlessly repeat coordinates.

compileBlocks :: Monad m => Enumeratee [Calls] [GenoCallBlock] m a
compileBlocks = convStream $ do
        c1 <- headStream
        tailBlock (p_refseq c1) (p_pos c1) (p_pos c1) . (:[]) $! pack c1
  where
    tailBlock !rs !p0 !po acc = do
        mc <- peekStream
        case mc of
            Just c1 | rs == p_refseq c1 && po+1 == p_pos c1 && po - p0 < 65536 -> do
                    _ <- headStream
                    tailBlock rs p0 (po+1) . (:acc) $! pack c1

            _ -> return [ GenoCallBlock
                    { reference_name = rs
                    , start_position = p0
                    , called_sites   = reverse acc } ]

    pack c1 = rlist indel_variants `seq` GenoCallSite{..}
      where
        Snp_GLs snp_pls !ref_allele = p_snp_pile c1

        !snp_stats         = p_snp_stat c1
        !indel_stats       = p_indel_stat c1
        !snp_likelihoods   = compact_likelihoods snp_pls
        !indel_likelihoods = compact_likelihoods $ fst $ p_indel_pile c1
        !indel_variants    = snd $ p_indel_pile c1

        rlist [] = ()
        rlist (x:xs) = x `seq` rlist xs


output_avro :: Handle -> Refs -> Iteratee [Calls] IO ()
output_avro hdl refs = compileBlocks =$
                       writeAvroContainer ContainerOpts{..} =$
                       mapChunksM_ (S.hPut hdl)
  where
    objects_per_block = 16
    filetype_label = "Genotype Likelihoods V0.1"
    initial_schemas = H.singleton "Refseq" $
        object [ "type" .= String "enum"
               , "name" .= String "Refseq"
               , "symbols" .= Array
                    (V.fromList . map (String . T.decodeUtf8 . sq_name) $ F.toList refs) ]
    meta_info = H.singleton "biohazard.refseq_length" $
                S.concat $ BL.toChunks $ encode $ Array $ V.fromList
                [ Number (fromIntegral (sq_length s)) | s <- F.toList refs ]


maxD :: Int
maxD = 64

-- | Parameter estimation for a single sample.  The parameters are
-- divergence and heterozygosity.  We tabulate the data here and do the
-- estimation afterwards.  Returns the product of the
-- parameter-independent parts of the likehoods and the histogram
-- indexed by D and H (see @genotyping.pdf@ for details).
tabulateSingle :: (Functor m, MonadIO m) => Iteratee [Calls] m DivTable
tabulateSingle = do
    tab <- liftIO $ M.replicate (12 * maxD * maxD) (0 :: Int)
    DivTable <$> foldStreamM (\acc -> accum tab acc . p_snp_pile) (0 :: Double)
             <*> liftIO (U.unsafeFreeze tab)
  where
    -- We need GL values for the invariant, the three homozygous variant
    -- and the three single-event heterozygous variant cases.  The
    -- ordering is like in BCF, with the reference first.
    -- Ref ~ A ==> PL ~ AA, AC, CC, AG, CG, GG, AT, CT, GT, TT
    {-# INLINE accum #-}
    accum !tab !acc (Snp_GLs !gls !ref)
        | U.length gls /= 10                   = error "Ten GL values expected for SNP!"      -- should not happen
        | ref `elem` [nucsC,nucsG]             = accum' 0 tab acc gls
        | ref `elem` [nucsA,nucsT]             = accum' 6 tab acc gls
        | otherwise                            = return acc                                   -- unknown reference

    -- The simple 2D table didn't work, it lacked resolution in some
    -- cases.  We make six separate tables instead so we can store two
    -- differences with good resolution in every case.
    {-# INLINE accum' #-}
    accum' refix !tab !acc !gls
        | g_RR >= g_RA && g_RA >= g_AA = store 0 g_RR g_RA g_AA
        | g_RR >= g_AA && g_AA >= g_RA = store 1 g_RR g_AA g_RA
        | g_RA >= g_RR && g_RR >= g_AA = store 2 g_RA g_RR g_AA
        | g_RA >= g_AA && g_AA >= g_RR = store 3 g_RA g_AA g_RR
        | g_RR >= g_RA                 = store 4 g_AA g_RR g_RA
        | otherwise                    = store 5 g_AA g_RA g_RR

      where
        g_RR = unPr $  U.unsafeIndex gls 0
        g_RA = unPr $ (U.unsafeIndex gls 1 + U.unsafeIndex gls 3 + U.unsafeIndex gls 6) / 3
        g_AA = unPr $ (U.unsafeIndex gls 2 + U.unsafeIndex gls 5 + U.unsafeIndex gls 9) / 3

        store t a b c = do let d1 = min (maxD-1) . round $ a - b
                               d2 = min (maxD-1) . round $ b - c
                               ix = (t + refix) * maxD * maxD + d1 * maxD + d2
                           liftIO $ M.read tab ix >>= M.write tab ix . succ
                           return $! acc + a

