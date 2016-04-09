{-# LANGUAGE BangPatterns, RecordWildCards, OverloadedStrings, FlexibleContexts #-}
-- Genotype calling for a single individual.  The only parameters needed
-- are the (prior) probabilities for a heterozygous or homozygous variant.
--
-- (This can be extended easily into a caller for a homogenous
-- population where individuals are assumed to be randomly related (i.e.
-- not family).  In this case, the prior is the allele frequency
-- spectrum, the call would be the set(!) of genotypes that has maximum
-- posterior probability.  Computation is possible in quadratic time and
-- linear space using a DP scheme; see Heng Li's paper for details.)
--
-- What's the output format?  Fasta or Fastq could be useful in limited
-- circumstances, else BCF (not VCF) would be canonical.  Or maybe BCF
-- restricted to variant sites.  Or BCF restricted to sites not known to
-- be reference.  People will come up with filters for sure...

import Bio.Base
import Bio.Bam
import Bio.Bam.Pileup
import Bio.Genocall.AvroFile
import Bio.Genocall.Metadata
import Bio.Iteratee.Builder
import Control.Applicative
import Control.Exception                ( bracket )
import Control.Monad
import Data.Avro
import Data.Bits
import Data.Foldable                    ( toList, foldMap )
import Data.MiniFloat
import Data.String
import Data.Text                        ( Text, unpack )
import Foreign.Ptr                      ( castPtr )
import System.Console.GetOpt
import System.Environment
import System.Exit
import System.FilePath
import System.Posix.IO

import qualified Data.ByteString.Char8          as S
import qualified Data.ByteString.Unsafe         as S
import qualified Data.HashMap.Strict            as H
import qualified Data.Vector.Unboxed            as U
import qualified System.IO                      as IO

data Conf = Conf {
    conf_metadata    :: FilePath,
    conf_ploidy      :: S.ByteString -> Int,
    conf_report      :: String -> IO () }

defaultConf :: Conf
defaultConf = Conf (error "no metadata file specified") (const 2) (\_ -> return ())

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "c"  ["config"]              (ReqArg set_conf "FILE") "Set name of json config file to FILE",
    Option "1"  ["haploid-chromosomes"] (ReqArg   set_hap "PRF") "Targets starting with PRF are haploid",
    Option "2"  ["diploid-chromosomes"] (ReqArg   set_dip "PRF") "Targets starting with PRF are diploid",
    Option "v"  ["verbose"]             (NoArg       be_verbose) "Print more diagnostics",
    Option "h?" ["help","usage"]        (NoArg       disp_usage) "Display this message" ]
  where
    disp_usage _ = do pn <- getProgName
                      let blah = "Usage: " ++ pn ++ " [OPTION...] [AVRO-FILE]"
                      putStrLn $ usageInfo blah options
                      exitFailure

    be_verbose       c = return $ c { conf_report = IO.hPutStrLn stderr }
    set_conf      fn c = return $ c { conf_metadata = fn }

    set_hap a c = return $ c { conf_ploidy = \chr -> if S.pack a `S.isPrefixOf` chr then 1 else conf_ploidy c chr }
    set_dip a c = return $ c { conf_ploidy = \chr -> if S.pack a `S.isPrefixOf` chr then 2 else conf_ploidy c chr }


main :: IO ()
main = do
    (opts, samples, errs) <- getOpt Permute options <$> getArgs
    unless (null errs) $ mapM_ (IO.hPutStrLn stderr) errs >> exitFailure
    when (null samples) $ IO.hPutStrLn stderr "need (at least) one sample name" >> exitFailure
    conf <- foldl (>>=) (return defaultConf) opts
    forM_ samples $ \sample -> do
        meta <- readMetadata (fromString (conf_metadata conf))

        case H.lookup (fromString sample) meta of
            Nothing -> IO.hPutStrLn stderr $ "unknown sample " ++ show sample
            Just smp -> main' conf sample smp


main' :: Conf -> String -> Sample -> IO ()
main' Conf{..} sample_name Sample{..} =
    case sample_divergences of
        Just (prior_div : prior_het : _prior_het2 : more) ->
            let prior_indel = case more of [] -> prior_div * 0.1 ; p : _ -> p in
            bracket (openFd (mkpath sample_bcf_file) WriteOnly (Just 0o666) defaultFileFlags) closeFd $ \ofd ->
            enumFile defaultBufSize (mkpath sample_avro_file) >=> run                                 $
            joinI $ readAvroContainer                                                                 $ \av_meta ->
            joinI $ progressPos getpos "calling at " conf_report (getRefseqs av_meta)                 $
            bcf_to_fd ofd (getRefseqs av_meta) [fromString sample_name]
                          (call (prior_div/3) prior_het, call prior_indel prior_het)

        _ -> fail $ show sample_name ++ " is missing divergence information"
  where
    call :: Double -> Double -> U.Vector (Prob' Float) -> Int
    call prior prior_h lks = U.maxIndex . U.zipWith (*) lks $
                             U.replicate (U.length lks) (toProb . realToFrac $ prior_h * prior)
                             U.// [ (0, toProb . realToFrac $ 1-prior) ]
                             U.// [ (i, toProb . realToFrac $ (1-prior_h) * prior)
                                  | i <- takeWhile (< U.length lks) (scanl (+) 2 [3..]) ]

    mkpath :: Text -> FilePath
    mkpath t = takeDirectory conf_metadata </> unpack t

    getpos :: GenoCallBlock -> (Refseq, Int)
    getpos b = (reference_name b, start_position b)


type CallFuncs = (U.Vector (Prob' Float) -> Int, U.Vector (Prob' Float) -> Int)

-- | Generates BCF and writes it to a 'Handle'.  For the necessary VCF
-- header, we get the /names/ of the reference sequences from the Avro
-- schema and the /lengths/ from the biohazard.refseq_length entry in
-- the meta data.
bcf_to_fd :: MonadIO m => Fd -> Refs -> [S.ByteString] -> CallFuncs -> Iteratee [GenoCallBlock] m ()
bcf_to_fd hdl refs smps call_fns = toBcf refs smps call_fns ><> encodeBgzfWith 9 =$ mapChunksM_ (liftIO . put)
  where
    put s = S.unsafeUseAsCStringLen s $ \(p,l) ->
            fdWriteBuf hdl (castPtr p) (fromIntegral l)

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


-- XXX Ploidy is being ignored.
toBcf :: Monad m => Refs -> [S.ByteString] -> CallFuncs -> Enumeratee [GenoCallBlock] Push m r
toBcf refs smps (snp_call, indel_call) = eneeCheckIfDone go
  where
    go  k = mapChunks (foldMap encode) . k $ Chunk hdr

    hdr   = pushByteString "BCF\2\2" <> setMark <>
            vcf_header refs smps <> pushByte 0 <> endRecord

    encode :: GenoCallBlock -> Push
    encode GenoCallBlock{..} = mconcat $ zipWith (encode1 reference_name) [start_position..] called_sites

    encode1 ref pos site =
        encodesnp ref pos site <>
        case indel_variants site of [ ] -> mempty
                                    [_] -> mempty
                                    _   -> encodeindel ref pos site

    encodesnp ref pos site = encode' ref pos (map S.singleton alleles) site (snp_stats site) snp_call
      where
        -- Permuting the reference allele to the front sucks.  Since
        -- there are only four possibilities, I'm not going to bother
        -- with an algorithm and just open-code it.
        alleles | ref_allele site == nucsT = "TACG"
                | ref_allele site == nucsG = "GACT"
                | ref_allele site == nucsC = "CAGT"
                | otherwise                = "ACGT"

    encodeindel ref pos site = encode' ref pos alleles site (indel_stats site) indel_call
      where
        -- We're looking at the indel /after/ the current position.
        -- That's sweet, because we can just prepend the current
        -- reference base and make bcftools and friends happy.  Longest
        -- reported deletion becomes the reference allele.  Others may
        -- need padding.
        rallele = snd $ maximum [ (U.length r, r) | IndelVariant (V_Nucs r) _ <- indel_variants site ]
        alleles = [ S.pack $ show (ref_allele site) ++ show (U.toList a) ++ show (U.toList $ U.drop (U.length r) rallele)
                  | IndelVariant (V_Nucs r) (V_Nuc a) <- indel_variants site ]

    encode' ref pos alleles GenoCallSite{..} CallStats{..} do_call =
        setMark <> setMark <>           -- remember space for two marks
        b_share <> endRecordPart1 <>    -- store 1st length and 2nd mark
        b_indiv <> endRecordPart2       -- store 2nd length
      where
        b_share = pushWord32 (unRefseq ref) <>
                  pushWord32 (fromIntegral pos) <>
                  pushWord32 0 <>                                   -- rlen?!  WTF?!
                  pushFloat gq <>                                   -- QUAL
                  pushWord16 2 <>                                   -- ninfo
                  pushWord16 (fromIntegral $ length alleles) <>     -- n_allele
                  pushWord32 0x04000001 <>                          -- n_fmt, n_sample
                  pushByte 0x07 <>                                  -- variant name (empty)
                  foldMap typed_string alleles <>                   -- alleles
                  pushByte 0x01 <>                                  -- FILTER (an empty vector)

                  pushByte 0x11 <> pushByte 0x01 <>                 -- INFO key 0 (MQ)
                  pushByte 0x11 <> pushByte rms_mapq <>             -- MQ, typed word8
                  pushByte 0x11 <> pushByte 0x02 <>                 -- INFO key 0 (MQ0)
                  pushByte 0x12 <> pushWord16 (fromIntegral reads_mapq0)   -- MQ0

        b_indiv = pushByte 0x01 <> pushByte 0x03 <>                 -- FORMAT key GT
                  pushByte 0x21 <>                                  -- two uint8s for GT
                  pushByte (2 + 2 * fromIntegral g) <>              -- actual GT
                  pushByte (2 + 2 * fromIntegral h) <>

                  pushByte 0x01 <> pushByte 0x04 <>                 -- FORMAT key DP
                  pushByte 0x12 <>                                  -- one uint16 for DP
                  pushWord16 (fromIntegral read_depth) <>           -- depth

                  pushByte 0x01 <> pushByte 0x05 <>                 -- FORMAT key PL
                  ( let l = U.length lks in if l < 15
                    then pushByte (fromIntegral l `shiftL` 4 .|. 2)
                    else pushWord16 0x03F2 <> pushWord16 (fromIntegral l) ) <>
                  pl_vals <>                                        -- vector of uint16s for PLs

                  pushByte 0x01 <> pushByte 0x06 <>                 -- FORMAT key GQ
                  pushByte 0x11 <> pushByte gq'                     -- uint8, genotype

        rms_mapq = round $ sqrt (fromIntegral sum_mapq_squared / fromIntegral read_depth :: Double)
        typed_string s | S.length s < 15 = pushByte (fromIntegral $ (S.length s `shiftL` 4) .|. 0x7) <> pushByteString s
                       | otherwise       = pushByte 0xF7 <> pushByte 0x03 <> pushWord32 (fromIntegral $ S.length s) <> pushByteString s

        pl_vals = U.foldr ((<>) . pushWord16 . round . max 0 . min 0x7fff . (*) (-10/log 10) . unPr . (/ lks U.! maxidx)) mempty lks

        lks = U.map (Pr . negate . mini2float) snp_likelihoods :: U.Vector (Prob' Float)
        maxidx = U.maxIndex lks

        gq = -10 * unPr (U.sum (U.ifilter (\i _ -> i /= maxidx) lks) / U.sum lks) / log 10
        gq' = round . max 0 . min 127 $ gq

        callidx = do_call lks
        h = length $ takeWhile (<= callidx) $ scanl (+) 1 [2..]
        g = callidx - h * (h+1) `div` 2



