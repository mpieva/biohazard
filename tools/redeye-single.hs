-- Genotype calling for a single individual on a subset of the genome.
-- It takes the appropriate set of .av files as input, along with priors
-- for divergence, indel rate, etc.  These also have sensible defaults.
-- Output is a BCF file.
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

import Bio.Bam
import Bio.Bam.Pileup
import Bio.Genocall.AvroFile
import Bio.Iteratee.Builder
import Bio.Prelude
import Bio.Util.Regex                   ( regComp, regMatch )
import Data.Avro
import Data.MiniFloat
import System.Console.GetOpt
import System.FilePath           hiding ( combine )
import System.Random

import qualified Data.ByteString.Char8          as S
import qualified Data.Vector.Unboxed            as U
import qualified System.IO                      as IO

type Reporter = String -> IO ()

-- XXX  conf_ploidy is completely ignored.  Should it be?
data Conf = Conf {
    conf_report     :: Reporter,
    conf_random     :: Maybe StdGen,
    conf_ploidy     :: String -> Int,
    sample_name     :: Conf -> String,
    outfile         :: FilePath,
    prior_div       :: Double,
    prior_het       :: Conf -> Double,
    prior_indel     :: Conf -> Double,
    prior_het_indel :: Conf -> Double }


defaultConf :: Conf
defaultConf = Conf { conf_report = \_ -> return ()
                   , conf_random = Nothing
                   , conf_ploidy = const 2
                   , sample_name = takeBaseName . outfile
                   , outfile = "-"
                   , prior_div = 0.001
                   , prior_het = (*) 2 . prior_div
                   , prior_indel = (*) 0.1 . prior_div
                   , prior_het_indel = \c -> prior_indel c c * prior_het c c / prior_div c }


options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o"  ["output"]            (ReqArg set_output "FILE") "Set output file pattern to FILE",
    Option "1"  ["haploid-chromosomes"] (ReqArg set_hap "REGEX") "Targets matching REGEX are haploid",
    Option "2"  ["diploid-chromosomes"] (ReqArg set_dip "REGEX") "Targets matching REGEX are diploid",
    Option "s"  ["sample-genotypes"]     (OptArg set_rnd "SEED") "Sample genotypes from posterior",
    Option "N"  ["name"]                (ReqArg set_name "NAME") "Set sample name to NAME",
    Option "d"  ["divergence"]           (ReqArg set_div "PROB") "Set probability of a hom. SNP to PROB",
    Option "D"  ["heterozygosity"]       (ReqArg set_het "PROB") "Set probability of a het. SNP to PROB",
    Option "i"  ["indel"]              (ReqArg set_indel "PROB") "Set probability of a hom. InDel to PROB",
    Option "I"  ["het-indel"]         (ReqArg set_hindel "PROB") "Set probability of a het. InDel to PROB",
    Option "v"  ["verbose"]             (NoArg       be_verbose) "Print more diagnostics",
    Option "h?" ["help","usage"]        (NoArg       disp_usage) "Display this message" ]
  where
    disp_usage _ = do pn <- getProgName
                      let blah = "Usage: " ++ pn ++ " [OPTION...] [AV-FILE ...]"
                      putStrLn $ usageInfo blah options
                      exitFailure

    be_verbose       c = return $ c { conf_report = IO.hPutStrLn stderr }

    set_hap        a c = return $ c { conf_ploidy = \ch -> if regMatch (regComp a) ch then 1 else conf_ploidy c ch }
    set_dip        a c = return $ c { conf_ploidy = \ch -> if regMatch (regComp a) ch then 2 else conf_ploidy c ch }

    set_output    fn c = return $ c { outfile     =       fn }
    set_name      nm c = return $ c { sample_name = const nm }

    set_rnd  Nothing c = newStdGen >>= \g -> return $ c { conf_random = Just g }
    set_rnd (Just a) c = readIO  a >>= \s -> return $ c { conf_random = Just (mkStdGen s) }

    set_div        a c = readIO a >>= \x -> return $ c { prior_div       =       x }
    set_het        a c = readIO a >>= \x -> return $ c { prior_het       = const x }
    set_indel      a c = readIO a >>= \x -> return $ c { prior_indel     = const x }
    set_hindel     a c = readIO a >>= \x -> return $ c { prior_het_indel = const x }


main :: IO ()
main = do
    (opts, infiles, errs) <- getOpt Permute options <$> getArgs
    unless (null errs) $ mapM_ (IO.hPutStrLn stderr) errs >> exitFailure
    conf@Conf{..} <- foldl (>>=) (return defaultConf) opts

    withOutputFd outfile                                                                                    $ \ofd ->
      concatAvs infiles >=> run                                                                             $ \av_meta ->
        progressPos (reference_name &&& start_position) "calling at " conf_report (getRefseqs av_meta)     =$
        bcf_to_fd ofd (getRefseqs av_meta) [fromString $ sample_name conf] conf_random
                      ( call (prior_div/3) (prior_het conf)
                      , call (prior_indel conf) (prior_het_indel conf) )


withOutputFd :: FilePath -> (Fd -> IO a) -> IO a
withOutputFd "-" k = k stdOutput
withOutputFd  f  k = do
    r <- withFd (f++".#") WriteOnly (Just 0o666) defaultFileFlags k
    rename (f++".#") f
    return r


concatAvs :: (MonadIO m, MonadMask m, Avro a) => [FilePath] -> Enumerator' AvroMeta [a] m b
concatAvs [        ] = \k -> enumHandle defaultBufSize stdin (readAvroContainer k) >>= run
concatAvs (fp0:fps0) = \k -> enum1 fp0 k >>= go fps0
  where
    enum1 "-" k1 = enumHandle defaultBufSize stdin (readAvroContainer k1) >>= run
    enum1  fp k1 = enumFile   defaultBufSize    fp (readAvroContainer k1) >>= run

    go [       ] = return
    go (fp1:fps) = enum1 fp1 . const >=> go fps


call :: Double -> Double -> U.Vector (Prob' Float) -> Maybe StdGen -> (Int,Maybe StdGen)
call prior_d prior_h lks gen = case gen of
    Nothing -> ( U.maxIndex ps, Nothing )
    Just  g -> (            ix, Just g' )
      where
        (p,g') = randomR (0, 1) g
        ix     = U.length $ U.takeWhile (<p) $ U.map fromProb $
                 U.postscanl (+) 0 $ U.map (/ U.sum ps) ps
  where
    ps = U.zipWith (*) lks $ U.replicate (U.length lks) (toProb . realToFrac $ prior_h * (1-prior_d))
                             U.// [ (0, toProb . realToFrac $ (1-prior_d) * (1-prior_h)) ]
                             U.// [ (i, toProb . realToFrac $ (1-prior_d) * prior_h)
                                  | i <- takeWhile (< U.length lks) (scanl (+) 2 [3..]) ]


-- | Generates BCF and writes it to a 'Handle'.  For the necessary VCF
-- header, we get the /names/ of the reference sequences from the Avro
-- schema and the /lengths/ from the biohazard.refseq_length entry in
-- the meta data.
bcf_to_fd :: MonadIO m => Fd -> Refs -> [S.ByteString] -> gen -> CallFuncs gen -> Iteratee [GenoCallBlock] m ()
bcf_to_fd hdl refs name gen callz =
    toBcf refs name callz gen ><> encodeBgzfWith 9 =$
    mapChunksM_ (liftIO . fdPut hdl)


-- A function from likelihoods to called index.  It's allowed to require
-- a random number generator.
type CallFunc gen = U.Vector (Prob' Float) -> gen -> (Int, gen)
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


-- XXX Ploidy is being ignored.
toBcf :: Monad m => Refs -> [S.ByteString] -> CallFuncs gen -> gen -> Enumeratee [GenoCallBlock] Push m r
toBcf refs smps (snp_call, indel_call) gen0 = eneeCheckIfDone go
  where
    go  k = eneeCheckIfDone (go2 gen0) . k $ Chunk hdr
    go2 g k = tryHead >>= \mblock -> case mblock of
                    Nothing -> return $ liftI k
                    Just bl -> let (p,g1) = encode bl g
                               in eneeCheckIfDone (go2 g1) . k $ Chunk p

    hdr   = pushByteString "BCF\2\2" <> setMark <>
            vcf_header refs smps <> pushByte 0 <> endRecord

    -- encode :: GenoCallBlock -> gen -> (Push, gen)
    encode GenoCallBlock{..} = combine start_position called_sites
      where
        combine !_ [    ] gen = (mempty, gen)
        combine !p (s:ss) gen = (p1 <> p2, g2)
          where
            (p1, g1) = encode1 reference_name p s gen
            (p2, g2) = combine (succ p) ss g1

    -- encode1 :: Refseq -> Int -> GenoCallSite -> gen -> (Push, gen)
    encode1 ref pos site g0 = (p1 <> p2, g2)
      where
        (p1,g1) = encodeSNP site ref pos snp_call g0
        (p2,g2) = case indel_variants site of
                    [ ] -> (mempty,g1)
                    [_] -> (mempty,g1)
                    v:_ | U.null d && U.null i -> encodeIndel site ref pos indel_call g1
                        | otherwise            -> error "First indel variant should always be the reference."
                      where
                        IndelVariant (V_Nucs d) (V_Nuc i) = v


encodeSNP :: GenoCallSite -> Refseq -> Int -> CallFunc gen -> gen -> (Push, gen)
encodeSNP site = encodeVar (map S.singleton alleles) (snp_likelihoods site) (snp_stats site)
  where
    -- Permuting the reference allele to the front sucks.  Since
    -- there are only four possibilities, I'm not going to bother
    -- with an algorithm and just open-code it.
    alleles | ref_allele site == nucsT = "TACG"
            | ref_allele site == nucsG = "GACT"
            | ref_allele site == nucsC = "CAGT"
            | otherwise                = "ACGT"

encodeIndel :: GenoCallSite -> Refseq -> Int -> CallFunc gen -> gen -> (Push, gen)
encodeIndel site = encodeVar alleles (indel_likelihoods site) (indel_stats site)
  where
    -- We're looking at the indel /after/ the current position.
    -- That's sweet, because we can just prepend the current
    -- reference base and make bcftools and friends happy.  Longest
    -- reported deletion becomes the reference allele.  Others may
    -- need padding.
    rallele = snd $ maximum [ (U.length r, r) | IndelVariant (V_Nucs r) _ <- indel_variants site ]
    alleles = [ S.pack $ showNucleotides (ref_allele site) : show (U.toList a) ++ show (U.toList $ U.drop (U.length r) rallele)
              | IndelVariant (V_Nucs r) (V_Nuc a) <- indel_variants site ]

encodeVar :: [S.ByteString] -> U.Vector Mini -> CallStats -> Refseq -> Int -> CallFunc gen -> gen -> (Push, gen)
encodeVar alleles likelihoods CallStats{..} ref pos do_call gen =
    ( setMark <> setMark <>           -- remember space for two marks
      b_share <> endRecordPart1 <>    -- store 1st length and 2nd mark
      b_indiv <> endRecordPart2       -- store 2nd length
    , gen' )
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

    lks = U.map (Pr . negate . mini2float) likelihoods :: U.Vector (Prob' Float)
    maxidx = U.maxIndex lks

    gq = -10 * unPr (U.sum (U.ifilter (\i _ -> i /= maxidx) lks) / U.sum lks) / log 10
    gq' = round . max 0 . min 127 $ gq

    (callidx, gen') = do_call lks gen
    h = length $ takeWhile (<= callidx) $ scanl (+) 1 [2..]
    g = callidx - h * (h+1) `div` 2


