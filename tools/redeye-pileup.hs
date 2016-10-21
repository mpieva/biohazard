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
import Bio.Genocall.Estimators
import Bio.Prelude
import Bio.Util.Pretty
import Data.Binary.Serialise.CBOR.SequenceFiles
import System.Console.GetOpt
import System.IO

import qualified Data.Binary                    as Bin
import qualified Data.ByteString.Lazy           as BL
import qualified Data.Vector.Storable           as VS
import qualified Data.Vector.Unboxed            as U
import qualified Data.Vector.Unboxed.Mutable    as M
import qualified Data.Sequence                  as Z

data Conf = Conf {
    conf_output :: FilePath,
    conf_libs   :: [Lib],
    conf_chrom  :: String,
    conf_theta  :: Maybe Double,
    conf_report :: Reporter,
    conf_table  :: Maybe FilePath }

-- | We pair libraries with their appropriate damage models.  Right now,
-- we specify the damage model (\"-D\", complete data structure in
-- 'gparse' format), then list the libraries it applies to.
defaultConf :: Conf
defaultConf = Conf { conf_output = error "no output file"
                   , conf_libs   = []
                   , conf_chrom  = ""
                   , conf_theta  = Nothing
                   , conf_report = \_ -> return ()
                   , conf_table  = Nothing }

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o"  ["output"] (ReqArg set_output "FILE") "Set output filename to FILE",
    Option "c"  chrom      (ReqArg set_chrom  "NAME") "Restrict to chromosome NAME",
    Option "D"  ["damage"] (ReqArg set_dmg   "PARMS") "Set damage parameters to PARMS",
    Option "t"  dep_param  (ReqArg set_theta  "FRAC") "Set dependency coefficient to FRAC (\"N\" to turn off)",
    Option "T"  ["table"]  (ReqArg want_table "FILE") "Print table for divergence estimation to FILE",
    Option "v"  ["verbose"]        (NoArg be_verbose) "Print more diagnostics",
    Option "h?" ["help","usage"]   (NoArg disp_usage) "Display this message" ]
  where
    dep_param = ["theta","dependency-coefficient"]
    chrom     = ["chromosome","region"]

    disp_usage    _ = do pn <- getProgName
                         let blah = "Usage: " ++ pn ++ " [[OPTION...] [BAM-FILE...] ...]"
                         putStrLn $ usageInfo blah options
                         exitSuccess

    be_verbose    c = return $ c { conf_report = hPutStrLn stderr }
    want_table fp c = return $ c { conf_table  = Just fp }

    set_theta "N" c = return $ c { conf_theta = Nothing }
    set_theta   a c = (\t  ->  c { conf_theta = Just  t }) <$> readIO a

    set_dmg     a c = let upd d = return $ c { conf_libs = Lib [] (NewDamage d) : conf_libs c }
                      in either fail upd . pparse . fromString $ a

    set_chrom   a c = return $ c { conf_chrom  = a }
    set_output fn c = return $ c { conf_output = fn }


main :: IO ()
main = do
    (opts, [], errs) <- let put_file  f c = return $ c { conf_libs = put_file1 f $ conf_libs c }
                            put_file1 f (Lib fs d : ls) =   Lib (f:fs) d : ls
                            put_file1 f [             ] = [ Lib [f] UnknownDamage ]
                        in getOpt (ReturnInOrder put_file) options <$> getArgs

    Conf{..} <- foldl (>>=) (return defaultConf) opts
    unless (null errs) $ mapM_ (hPutStrLn stderr) errs >> exitFailure

    (tab,()) <- withFile (conf_output ++ ".#") WriteMode                                        $ \ohdl ->
                mergeLibraries conf_report (reverse conf_libs) conf_chrom >=> run               $ \hdr ->
                progressPos (\(rs, p, _) -> (rs, p)) "GT call at " conf_report (meta_refs hdr) =$
                pileup                                                                         =$
                mapStream (calls conf_theta)                                                   =$
                zipStreams tabulateSingle (output_cbor ohdl $ meta_refs hdr)

    rename (conf_output ++ ".#") conf_output

    maybe (return ()) (BL.hPut stdout . Bin.encode) conf_table

    (de1,de2) <- estimateSingle tab
    putStrLn $ unlines $
            showRes (point_est de1) :
            [ "[ " ++ showRes u ++ " .. " ++ showRes v ++ " ]" | (u,v) <- conf_region de1 ] ++
            [] : showRes (point_est de2) :
            [ "[ " ++ showRes u ++ " .. " ++ showRes v ++ " ]" | (u,v) <- conf_region de2 ]
  where
    showRes     [dv,h] = "D  = " ++ showFFloat (Just 6) dv ", " ++
                         "H  = " ++ showFFloat (Just 6) h ""
    showRes [dv,hs,hw] = "D  = " ++ showFFloat (Just 6) dv ", " ++
                         "Hs = " ++ showFFloat (Just 6) hs ", " ++
                         "Hw = " ++ showFFloat (Just 6) hw ""
    showRes          _ = error "Wtf? (showRes)"


mergeLibraries :: (MonadIO m, MonadMask m)
               => Reporter -> [Lib] -> String
               -> Enumerator' BamMeta [PosPrimChunks] m b
mergeLibraries       _ [    ]    _ = error "Need at least one library."
mergeLibraries  report [ l  ] mrgn = enumLibrary report l mrgn
mergeLibraries  report (l:ls) mrgn = mergeEnums' (mergeLibraries report ls mrgn) (enumLibrary report l mrgn) mm
  where
    mm _ = mergeSortStreams $ \(rs1, p1, _) (rs2, p2, _) -> if (rs1, p1) < (rs2, p2) then Less else NotLess

data Lib = Lib [String] (GenDamageParameters U.Vector Double)
type Reporter = String -> IO ()

enumLibrary :: (MonadIO m, MonadMask m)
            => Reporter -> Lib -> String -> Enumerator' BamMeta [PosPrimChunks] m b
enumLibrary report (Lib fs mdp) mrgn output = do
    let (msg, dmg) = case mdp of UnknownDamage -> ("no damage model",                               noDamage)
                                 OldDamage  dp -> ("universal damage parameters " ++ show dp,  univDamage dp)
                                 NewDamage ndp -> ("empirical damage parameters " ++ show ndp, empDamage ndp)

    liftIO . report $ "using " ++ msg ++ " for " ++ show fs

    mergeInputRgns combineCoordinates mrgn (reverse fs)
        $== takeWhileE (isValidRefseq . b_rname . unpackBam)
        $== mapMaybeStream (\br ->
                let b = unpackBam br
                    m = dmg (isReversed b) (VS.length (b_qual b))
                in decompose (map m [0..]) br)
        $ output

mergeInputRgns :: (MonadIO m, MonadMask m)
               => (BamMeta -> Enumeratee [BamRaw] [BamRaw] (Iteratee [BamRaw] m) a)
               -> String -> [FilePath] -> Enumerator' BamMeta [BamRaw] m a
mergeInputRgns  _   _ [        ] = \k -> return (k mempty)
mergeInputRgns (?) ""      fps   = mergeInputs (?) fps
mergeInputRgns (?) rs (fp0:fps0) = go fp0 fps0
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

compileBlocks :: Monad m => Enumeratee [Calls] [GenoFileRec] m a
compileBlocks = unfoldConvStream conv_one (Refseq 0, 0)
  where
    conv_one (!rs,!po) = do
        cc <- headStream

        let Snp_GLs snp_pls ref_allele = p_snp_pile cc
            snp_stats         = p_snp_stat cc
            indel_stats       = p_indel_stat cc
            snp_likelihoods   = compact_likelihoods snp_pls
            indel_likelihoods = compact_likelihoods $ fst $ p_indel_pile cc
            indel_variants    = snd $ p_indel_pile cc

            blk_hdr = if rs == p_refseq cc && po == p_pos cc
                            then id
                            else (:) ( GenoFileBlock (GenoCallBlock
                                     { reference_name = p_refseq cc, start_position = p_pos cc }))
        return ( (p_refseq cc, p_pos cc + 1), blk_hdr [ GenoFileSite GenoCallSite{..} ] )


output_cbor :: Handle -> Refs -> Iteratee [Calls] IO ()
output_cbor hdl refs =
    ( lift . enumPure1Chunk [GenoFileHeader (GenoHeader 0 refs)]    >=>
      compileBlocks                                                 >=>
      lift . enumPure1Chunk [GenoFileFooter GenoFooter] )            =$
    mapChunksM_ (hPutBinaryFileSequence hdl)


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

