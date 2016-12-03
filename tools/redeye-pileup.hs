{-# LANGUAGE TypeSynonymInstances, FlexibleInstances, CPP #-}
-- Command line driver for simple genotype calling.  We have three
-- separate steps:  Pileup from a BAM file (or multiple merged files) to
-- produce likelihoods (and some auxillary statistics).  These are
-- written into an CBOR sequence file.  Next we need to estimate parameters,
-- in the simplest case divergence and heterozygosity.  We can save some
-- time by fusing this with the first step.  The final step is calling
-- bases by scnaning the CBOR sequence file and applying some model, and
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
--
-- XXX  In here, we rely on read groups beind declared sensibly.  We
-- will estimate one damage/error model per read group.  I'm not sure
-- doing it per library would be more sensible, but then we'd have to
-- rely on yet another set of annotations and/or headers being present.

import Bio.Adna
import Bio.Bam
import Bio.Bam.Pileup
import Bio.Genocall
import Bio.Genocall.Estimators
import Bio.Genocall.LkFile
import Bio.Prelude
import Data.Aeson
import Data.Binary.Serialise.CBOR.SequenceFiles
import Data.Text.Encoding                       ( encodeUtf8 )
import System.Console.GetOpt
import System.IO

import qualified Data.Binary                    as Bin
import qualified Data.ByteString                as BS
import qualified Data.ByteString.Lazy           as BL
import qualified Data.HashMap.Strict            as H
import qualified Data.Sequence                  as Z
import qualified Data.Vector                    as V

data Conf = Conf {
    conf_output :: FilePath,
    conf_dmg    :: HashMap Bytes SubstModel,
    conf_chrom  :: String,
    conf_theta  :: Maybe Double,
    conf_report :: String -> IO (),
    conf_table  :: Maybe FilePath }

#if !MIN_VERSION_aeson(1,0,0)
instance FromJSON a => FromJSON (HashMap Bytes a) where
    parseJSON = withObject "hashmap" $ \o ->
                H.fromList <$> mapM (\(k,v) -> (,) (encodeUtf8 k) <$> parseJSON v) (H.toList o)
#else
instance FromJSONKey BS.ByteString where
    fromJSONKey = FromJSONKeyText encodeUtf8
    fromJSONKeyList = undefined -- XXX whatever
#endif

instance FromJSON Mat44D
instance FromJSON SubstModel

-- | We map read groups to damage models.  The set of damage models is
-- supplied in a JSON file.
defaultConf :: Conf
defaultConf = Conf { conf_output = error "no output file"
                   , conf_dmg    = mempty
                   , conf_chrom  = ""
                   , conf_theta  = Nothing
                   , conf_report = \_ -> return ()
                   , conf_table  = Nothing }

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o"  ["output"] (ReqArg set_output "FILE") "Set output filename to FILE",
    Option "c"  chrom      (ReqArg set_chrom  "NAME") "Restrict to chromosome NAME",
    Option "D"  ["damage"] (ReqArg set_dmg    "FILE") "Read damage model from FILE",
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

    set_dmg     a c = BL.readFile a >>= \s -> case eitherDecode' s of
                            Left err -> error err
                            Right ds -> return $ c { conf_dmg = ds }

    set_chrom   a c = return $ c { conf_chrom  = a }
    set_output fn c = return $ c { conf_output = fn }


main :: IO ()
main = do
    (opts, files, errs) <- getOpt Permute options <$> getArgs

    Conf{..} <- foldl (>>=) (return defaultConf) opts
    unless (null errs) $ mapM_ (hPutStrLn stderr) errs >> exitFailure

    (tab,()) <- withFile (conf_output ++ ".#") WriteMode                                        $ \ohdl ->
                mergeInputRgns combineCoordinates conf_chrom files >=> run                      $ \hdr ->
                takeWhileE (isValidRefseq . b_rname . unpackBam)                               =$
                concatMapStream (decompose_dmg_from conf_dmg)                                  =$
                progressPos (\(a,b,_)->(a,b)) "GT call at" (meta_refs hdr) 0x4000 conf_report  =$
                pileup                                                                         =$
                mapStream (calls conf_theta)                                                   =$
                zipStreams tabulateSingle (output_cbor ohdl $ meta_refs hdr)

    rename (conf_output ++ ".#") conf_output
    forM_ conf_table $ flip BL.writeFile $ Bin.encode tab

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


{-# INLINE decompose_dmg_from #-}
decompose_dmg_from :: HashMap Bytes SubstModel -> BamRaw -> [PosPrimChunks Mat44D]
decompose_dmg_from hm raw =
    decompose (model (H.lookup (extAsString "RG" (unpackBam raw)) hm)) raw
  where
    model              Nothing  _ _ = scalarMat 1
    model (Just SubstModel{..}) i r
        | i >= 0 &&   i  <  V.length  left_substs_fwd && not r = V.unsafeIndex left_substs_fwd    i
        | i <  0 && (-i) >= V.length right_substs_fwd && not r = V.unsafeIndex right_substs_fwd (-i-1)
        | not r                                                = middle_substs_fwd
        | i >= 0 &&   i  <  V.length  left_substs_rev          = V.unsafeIndex left_substs_rev    i
        | i <  0 && (-i) >= V.length right_substs_rev          = V.unsafeIndex right_substs_rev (-i-1)
        | otherwise                                            = middle_substs_rev


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

calls :: Maybe Double -> Pile Mat44D -> Calls
calls Nothing pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !s = simple_snp_call   $ uncurry (++) $ p_snp_pile pile
    !i = simple_indel_call $ p_indel_pile pile
    -- XXX this should be a cmdline option, if we ever look at qualities again
    -- fq = min 1 . (*) 1.333 . fromQual
    -- fq = fromQual

calls (Just _theta) _pile = error "Sorry, maq_snp_call is broken right now." -- XXX
{- calls (Just theta) pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !i = simple_indel_call $ p_indel_pile pile

    -- This lumps the two strands together
    -- !s = maq_snp_call theta $ uncurry (++) $ p_snp_pile pile -- XXX

    -- This treats them separately
    !s | r == r'    = Snp_GLs (U.zipWith (*) x y) r     -- same ref base (normal case): multiply
       | r == nucsN = Snp_GLs y r'                      -- forward ref is N, use backward call
       | otherwise  = Snp_GLs x r                       -- else use forward call (even if this is incorrect,
      where                                             -- there is nothing else we can do here)
        Snp_GLs x r  = maq_snp_call theta $ fst $ p_snp_pile pile
        Snp_GLs y r' = maq_snp_call theta $ snd $ p_snp_pile pile -}


-- | Serialize the results from genotype calling in a sensible way.  We
-- write a CBOR sequence file, but we use two different records, so we
-- don't need to endlessly repeat coordinates.

compileBlocks :: Monad m => Enumeratee [Calls] [GenoFileRec] m a
compileBlocks = unfoldConvStream conv_one (Refseq 0, 0)
  where
    conv_one (!rs,!po) = do
        cc <- headStream

        let ref_allele        = snp_refbase $ p_snp_pile cc
            snp_stats         = p_snp_stat cc
            indel_stats       = p_indel_stat cc
            snp_likelihoods   = compact_likelihoods $ snp_gls $ p_snp_pile cc
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


