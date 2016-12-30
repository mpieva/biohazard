{-# LANGUAGE TypeSynonymInstances, FlexibleInstances, CPP #-}
-- Command line driver for simple genotype calling.  We have two
-- separate steps:  We estimate parameters on a subset of the input,
-- then we pile up.  Output from pileup is a BCF file containing
-- likelihoods, posterior probabilities, and genotype calls.
-- Alternatively we could write a Heffalump, which has only genotype
-- calls.  We also write a table on the side, which can be used to
-- estimate divergence and heterozygosity per chromosome or genome wide.
--
-- The likelihoods depend on damage parameters and an error model,
-- otherwise they are 'eternal'.  (For the time being, it's probably
-- wise to go with the naïve error model.)  Technically, they also
-- depend on ploidy, but since only diploid organisms are interesting
-- right now, we fix that to two.  We pay some overhead on the sex
-- chromosomes, but the simplification is worth it.
--
-- About damage parameters:  Everything converged on an empirical model
-- where damage and sequencing error are both modelled as one position
-- specific substitution matrix.  The damage model is specific to a read
-- group, so read group annotations must be present.
--
-- Calling is always diploid, for maximum flexibility.  We don't really
-- support higher ploidies, so the worst damage is that we output an
-- overhead of 150% useless likelihood values for the sex chromosomes
-- and maybe estimate heterozygosity where there is none.
--
-- (This can be extended easily into a caller for a homogenous
-- population where individuals are assumed to be randomly related (i.e.
-- not family).  In this case, the prior is the allele frequency
-- spectrum, the call would be the set(!) of genotypes that has maximum
-- posterior probability.  Computation is possible in quadratic time and
-- linear space using a DP scheme; see Heng Li's paper for details.)

import Bio.Adna
import Bio.Bam
import Bio.Bam.Pileup
import Bio.Genocall
import Bio.Genocall.Estimators
import Bio.Iteratee.Builder
import Bio.Prelude
import Data.Aeson
import GHC.Float                                ( double2Float )
import System.Console.GetOpt
import System.FilePath
import System.Random

import qualified Data.Binary                    as Bin
import qualified Data.ByteString.Char8          as S
import qualified Data.ByteString.Lazy           as BL
import qualified Data.HashMap.Strict            as H
import qualified Data.Sequence                  as Z
import qualified Data.Vector                    as V
import qualified Data.Vector.Unboxed            as U

type Reporter = String -> IO ()

data Conf = Conf {
    conf_output     :: FilePath,
    conf_dmg        :: ExtModel,
    conf_chrom      :: String,
    conf_report     :: Reporter,
    conf_table      :: Maybe FilePath,
    conf_random     :: Maybe StdGen,
    sample_name     :: Conf -> String }


-- | We map read groups to damage models.  The set of damage models is
-- supplied in a JSON file.
defaultConf :: Conf
defaultConf = Conf { conf_output = error "no output file"
                   , conf_dmg    = ExtModel (DivEst [0.001,0.0005] []) Nothing (SubstModels mempty)
                   , conf_chrom  = ""
                   , conf_report = \_ -> return ()
                   , conf_table  = Nothing
                   , conf_random = Nothing
                   , sample_name = takeWhile (/='.') . takeFileName . conf_output }


options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o"  ["output"]              (ReqArg set_output "FILE") "Set output filename to FILE",
    Option "c"  ["chromosome","region"] (ReqArg set_chrom  "NAME") "Restrict to chromosome NAME",
    Option "D"  ["damage"]              (ReqArg set_dmg    "FILE") "Read damage model from FILE",
    Option "T"  ["table"]               (ReqArg want_table "FILE") "Print table for divergence estimation to FILE",
    Option "s"  ["sample-genotypes"]    (OptArg set_rnd    "SEED") "Sample genotypes from posterior",
    Option "N"  ["name"]                (ReqArg set_name   "NAME") "Set sample name to NAME",
    Option "v"  ["verbose"]             (NoArg         be_verbose) "Print more diagnostics",
    Option "h?" ["help","usage"]        (NoArg         disp_usage) "Display this message" ]
  where
    disp_usage    _ = do pn <- getProgName
                         let blah = "Usage: " ++ pn ++ " [[OPTION...] [BAM-FILE...] ...]"
                         putStrLn $ usageInfo blah options
                         exitSuccess

    be_verbose       c = return $ c { conf_report = hPutStrLn stderr }
    want_table    fp c = return $ c { conf_table  = Just fp }

    set_dmg        a c = BL.readFile a >>= \s -> case eitherDecode' s of
                            Left err -> error err
                            Right ds -> return $ c { conf_dmg = ds }

    set_chrom      a c = return $ c { conf_chrom  = a }
    set_output    fn c = return $ c { conf_output = fn }
    set_name      nm c = return $ c { sample_name = const nm }

    set_rnd  Nothing c = newStdGen >>= \g -> return $ c { conf_random = Just g }
    set_rnd (Just a) c = readIO  a >>= \s -> return $ c { conf_random = Just (mkStdGen s) }


main :: IO ()
main = do
    (opts, files, errs) <- getOpt Permute options <$> getArgs
    conf@Conf{..} <- foldl (>>=) (return defaultConf) opts
    unless (null errs) $ mapM_ (hPutStrLn stderr) errs >> exitFailure

    let symtab = H.fromList $ zip [ k | (k,_) <- case damage conf_dmg of SubstModels ms -> H.toList ms ] (map DmgToken [0..])
        smodel = V.fromList       [ v | (_,v) <- case damage conf_dmg of SubstModels ms -> H.toList ms ]

    let prior_div:prior_het:_ = point_est $ population conf_dmg
        callz = ( call $ SinglePop prior_div prior_het
                , call $ SinglePop (0.1 * prior_div) (0.1 * prior_het) )

    (tab,()) <- withOutputFd conf_output                                                         $ \ofd ->
                mergeInputRgns combineCoordinates conf_chrom files >=> run                       $ \hdr ->
                takeWhileE (isValidRefseq . b_rname . unpackBam)                                =$
                concatMapStream (decompose_dmg_from symtab)                                     =$
                progressPos (\(a,b,_,_)->(a,b)) "GT call at" (meta_refs hdr) 0x4000 conf_report =$
                pileup                                                                          =$
                mapStream (simple_calls smodel)                                                 =$
                zipStreams tabulateSingle
                           (toBcf (meta_refs hdr) [fromString $ sample_name conf] callz conf_random     =$
                            mapChunksM_ (liftIO . fdPut ofd))

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


withOutputFd :: FilePath -> (Fd -> IO a) -> IO a
withOutputFd "-" k = k stdOutput
withOutputFd  f  k = do
    r <- withFd (f++".#") WriteOnly (Just 0o666) defaultFileFlags k
    rename (f++".#") f
    return r

simple_calls :: V.Vector SubstModel -> Pile -> Calls
simple_calls dmods pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !s = simple_snp_call   get_dmg $ p_snp_pile pile
    !i = simple_indel_call get_dmg $ p_indel_pile pile

    get_dmg :: DmgToken -> Int -> Bool -> Mat44D
    get_dmg (DmgToken dt) = case dmods V.!? dt of
        Nothing -> error "shouldn't happen"
        Just sm -> lookupSubstModel sm

{-# INLINE decompose_dmg_from #-}
decompose_dmg_from :: HashMap Bytes DmgToken -> BamRaw -> [PosPrimChunks]
decompose_dmg_from hm raw =
    let rg = extAsString "RG" (unpackBam raw)
    in case H.lookup rg hm of
            Just tk -> decompose tk raw
            Nothing -> error $ "no substitution model for " ++ unpack rg


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

call :: SinglePop -> U.Vector Prob -> Maybe StdGen -> (Int,Maybe StdGen)
call priors lks gen = case gen of
    Nothing -> ( U.maxIndex ps, Nothing )
    Just  g -> (            ix, Just g' )
      where
        (p,g') = randomR (0, 1) g
        ix     = U.length $ U.takeWhile (<p) $ U.map fromProb $
                 U.init $ U.postscanl (+) 0 $ U.map (/ U.sum ps) ps
  where
    ps = single_pop_posterior priors 0 lks



-- A function from likelihoods to called index.  It's allowed to require
-- a random number generator.
type CallFunc gen = U.Vector Prob -> gen -> (Int, gen)
type CallFuncs gen = (CallFunc gen, CallFunc gen)

vcf_header :: Refs -> [S.ByteString] -> Push
vcf_header refs smps = foldr (\a b -> pushByteString a <> pushByte 10 <> b) mempty $
    [ "##fileformat=VCFv4.2"
    , "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"RMS mapping quality\">"
    , "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Number of MAPQ==0 reads covering this record\">"
    , "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    , "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"read depth\">"
    , "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"genotype likelihoods in deciban\">"
    , "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"conditional genotype quality in deciban\">" ] ++
    [ S.concat [ "##contig=<ID=", sq_name s, ",length=", S.pack (show (sq_length s)), ">" ] | s <- toList refs ] ++
    [ S.intercalate "\t" $ "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" : smps ]


toBcf :: MonadIO m => Refs -> [S.ByteString] -> CallFuncs gen -> gen -> Enumeratee [Calls] Bytes m r
toBcf refs smps (snp_call, indel_call) gen0 = eneeCheckIfDone go ><> encodeBgzfWith 6
  where
    go    k = eneeCheckIfDone (go2 gen0) . k $ Chunk hdr
    go2 g k = tryHead >>= \zz -> case zz of
                    Nothing -> return $ liftI k
                    Just cs -> let (p,g1) = encode1 cs g
                               in eneeCheckIfDone (go2 g1) . k $ Chunk p

    hdr     = pushByteString "BCF\2\2" <> setMark <>
              vcf_header refs smps <> pushByte 0 <> endRecord

    encode1 cs g0 = (p1 <> p2, g2)
      where
        (p1,g1) = encodeSNP cs snp_call g0
        (p2,g2) = case snd $ p_indel_pile cs of
                    [ ] -> (mempty,g1)
                    [_] -> (mempty,g1)
                    v:_ | U.null d && U.null i -> encodeIndel cs indel_call g1
                        | otherwise            -> error "First indel variant should always be the reference."
                      where
                        IndelVariant (V_Nucs d) (V_Nuc i) = v


encodeSNP :: Calls -> CallFunc gen -> gen -> (Push, gen)
encodeSNP cs = encodeVar (map S.singleton alleles) (gls `U.backpermute` U.fromList perm)
                         (p_snp_stat cs) (p_refseq cs) (p_pos cs)
  where
    Snp_GLs gls ref_allele = p_snp_pile cs

    -- Permuting the reference allele to the front in alleles and PLs
    -- sucks.  Since there are only four possibilities, I'm not going to
    -- bother with an algorithm and just open-code it.
    (alleles, perm) | ref_allele == nucsT = ("TACG", [9,6,0,7,1,2,8,3,4,5])
                    | ref_allele == nucsG = ("GACT", [5,3,0,4,1,2,8,6,7,9])
                    | ref_allele == nucsC = ("CAGT", [2,1,0,4,3,5,7,6,8,9])
                    | otherwise           = ("ACGT", [0,1,2,3,4,5,6,7,8,9])

encodeIndel :: Calls -> CallFunc gen -> gen -> (Push, gen)
encodeIndel cs = encodeVar alleles (fst $ p_indel_pile cs)
                           (p_indel_stat cs) (p_refseq cs) (p_pos cs)
  where
    Snp_GLs _ ref_allele = p_snp_pile cs

    -- We're looking at the indel /after/ the current position.  That's
    -- sweet, because we can just prepend the current reference base and
    -- make bcftools and friends happy.  Longest reported deletion
    -- becomes the reference allele.  Others may need padding.
    rallele = snd $ maximum [ (U.length r, r) | IndelVariant (V_Nucs r) _ <- snd $ p_indel_pile cs ]
    alleles = [ S.pack $ showNucleotides ref_allele : show (U.toList a) ++ show (U.toList $ U.drop (U.length r) rallele)
              | IndelVariant (V_Nucs r) (V_Nuc a) <- snd $ p_indel_pile cs ]

encodeVar :: [S.ByteString] -> U.Vector Prob -> CallStats -> Refseq -> Int -> CallFunc gen -> gen -> (Push, gen)
encodeVar alleles lks CallStats{..} ref pos do_call gen =
    ( setMark <> setMark <>           -- remember space for two marks
      b_share <> endRecordPart1 <>    -- store 1st length and 2nd mark
      b_indiv <> endRecordPart2       -- store 2nd length
    , gen' )
  where
    b_share = pushWord32 (unRefseq ref) <>
              pushWord32 (fromIntegral pos) <>
              pushWord32 0 <>                                   -- rlen?!  WTF?!
              pushFloat (double2Float gq) <>                    -- QUAL
              pushWord16 2 <>                                   -- ninfo
              pushWord16 (fromIntegral $ length alleles) <>     -- n_allele
              pushWord32 0x04000001 <>                          -- n_fmt, n_sample
              pushByte 0x07 <>                                  -- variant name (empty)
              foldMap typed_string alleles <>                   -- alleles
              pushByte 0x01 <>                                  -- FILTER (an empty vector)

              pushByte 0x11 <> pushByte 0x01 <>                 -- INFO key 0 (MQ)
              pushByte 0x11 <> pushByte rms_mapq <>             -- MQ, typed word8
              pushByte 0x11 <> pushByte 0x02 <>                 -- INFO key 1 (MQ0)
              pushByte 0x12 <> pushWord16 (fromIntegral reads_mapq0) -- MQ0

    b_indiv = pushByte 0x01 <> pushByte 0x03 <>                 -- FORMAT key 2 (GT)
              pushByte 0x21 <>                                  -- two uint8s for GT
              pushByte (2 + 2 * fromIntegral g) <>              -- actual GT
              pushByte (2 + 2 * fromIntegral h) <>

              pushByte 0x01 <> pushByte 0x04 <>                 -- FORMAT key 3 (DP)
              pushByte 0x12 <>                                  -- one uint16 for DP
              pushWord16 (fromIntegral read_depth) <>           -- depth

              pushByte 0x01 <> pushByte 0x05 <>                 -- FORMAT key 4 (PL)
              ( let l = U.length lks in if l < 15
                then pushByte (fromIntegral l `shiftL` 4 .|. 2)
                else pushWord16 0x02F2 <> pushWord16 (fromIntegral l) ) <>
              pl_vals <>                                        -- vector of uint16s for PLs

              pushByte 0x01 <> pushByte 0x06 <>                 -- FORMAT key 5 (GQ)
              pushByte 0x11 <> pushByte gq'                     -- uint8, genotype

    rms_mapq = round $ sqrt (fromIntegral sum_mapq_squared / fromIntegral read_depth :: Double)
    typed_string s | S.length s < 15 = pushByte (fromIntegral $ S.length s `shiftL` 4 .|. 0x7) <> pushByteString s
                   | otherwise       = pushByte 0xF7 <> pushByte 0x03 <> pushWord32 (fromIntegral $ S.length s) <> pushByteString s

    pl_vals = U.foldr ((<>) . pushWord16 . round . max 0 . min 0x7fff . (*) (-10/log 10) . unPr . (/ lks U.! maxidx)) mempty lks
    maxidx = U.maxIndex lks

    gq = -10 * unPr (U.sum (U.ifilter (\i _ -> i /= maxidx) lks) / U.sum lks) / log 10
    gq' = round . max 0 . min 127 $ gq

    (callidx, gen') = do_call lks gen
    h = length $ takeWhile (<= callidx) $ scanl (+) 1 [2..]
    g = callidx - h * (h+1) `div` 2
