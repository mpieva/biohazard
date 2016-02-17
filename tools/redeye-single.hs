{-# LANGUAGE BangPatterns, RecordWildCards, OverloadedStrings, FlexibleContexts #-}
-- Genotype calling for a single individual.  The only parameters needed
-- are the (prior) probabilities for a heterozygous or homozygous variant.
--
-- (This can be extended easily into a caller for a homogenous
-- population where individuals are assumed to be randomly related (i.e.
-- not family).  In this case, the prior is the allele frequency
-- spectrum, the call would be the set(!) of genotypes that has maximum
-- posterior probability.  Computation is possible in quadratic time and
-- linear space usind a DP scheme; see Heng Li's paper for details.)
--
-- What's the output format?  Fasta or Fastq could be useful in limited
-- circumstances, else BCF (not VCF) would be canonical.  Or maybe BCF
-- restricted to variant sites.  Or BCF restricted to sites not known to
-- be reference.  People will come up with filters for sure...

import Bio.Base
import Bio.Bam
import Bio.Bam.Pileup
import Bio.Genocall.AvroFile
import Bio.Iteratee.Builder
import Control.Applicative
import Control.Exception ( bracket )
import Control.Monad
import Data.Avro
import Data.Bits
import Data.Foldable ( toList, foldMap )
import Data.List ( scanl )
import Data.MiniFloat
import Data.Monoid
import Foreign.Ptr ( castPtr )
import System.Console.GetOpt
import System.Environment
import System.Exit
import System.Posix.IO

import qualified Data.ByteString.Char8          as S
import qualified Data.ByteString.Unsafe         as S
import qualified Data.Vector.Unboxed            as U
import qualified System.IO                      as IO

data Conf = Conf {
    conf_output      :: Output (),
    conf_sample      :: S.ByteString,
    conf_ploidy      :: S.ByteString -> Int,
    conf_report      :: String -> IO (),
    conf_prior_div   :: Prob' Float,
    conf_prior_het   :: Prob' Float,
    conf_prior_indel :: Prob' Float }

defaultConf :: Conf
defaultConf = Conf ($ bcf_to_fd stdOutput) "John_Doe" (const 2) (\_ -> return ())
                   (qualToProb $ Q 30) (2/3) (qualToProb $ Q 45)

options :: [OptDescr (Conf -> IO Conf)]
options = [ -- Maybe add FastQ output?  Or FastA?  Or padded FastA?  Or Nick FastA?
    Option "o" ["output", "bcf-output"]             (ReqArg set_bcf_out "FILE") "Write BCF output to FILE",
    Option "N" ["name","sample-name"]               (ReqArg  set_sample "NAME") "Set sample name to NAME",
    Option "1" ["haploid-chromosomes"]              (ReqArg     set_hap  "PRF") "Targets starting with PRF are haploid",
    Option "2" ["diploid-chromosomes"]              (ReqArg     set_dip  "PRF") "Targets starting with PRF are diploid",
    Option "D" ["prior-variant", "divergence"]      (ReqArg     set_div "PROB") "Set prior for a SNP variant to PROB",
    Option "I" ["prior-indel", "indel-rate"]        (ReqArg   set_indel "PROB") "Set prior for an indel variant to PROB",
    Option "H" ["prior-heterozygous", "heterozygosity"] (ReqArg set_het "PROB") "Set prior for a heterozygous variant to PROB",
    Option "v" ["verbose"]                          (NoArg          be_verbose) "Print more diagnostics",
    Option "h?" ["help","usage"]                    (NoArg          disp_usage) "Display this message" ]
  where
    disp_usage _ = do pn <- getProgName
                      let blah = "Usage: " ++ pn ++ " [OPTION...] [AVRO-FILE]"
                      putStrLn $ usageInfo blah options
                      exitFailure

    be_verbose       c = return $ c { conf_report = IO.hPutStrLn stderr }
    set_sample     a c = return $ c { conf_sample = S.pack a }
    set_bcf_out  "-" c = return $ c { conf_output = ($ bcf_to_fd stdOutput) }
    set_bcf_out   fn c = return $ c { conf_output = \k -> withFd fn WriteOnly (k . bcf_to_fd) }

    set_hap a c = return $ c { conf_ploidy = \chr -> if S.pack a `S.isPrefixOf` chr then 1 else conf_ploidy c chr }
    set_dip a c = return $ c { conf_ploidy = \chr -> if S.pack a `S.isPrefixOf` chr then 2 else conf_ploidy c chr }

    set_div        a c = (\p -> c { conf_prior_div   = toProb p }) <$> readIO a
    set_indel      a c = (\p -> c { conf_prior_indel = toProb p }) <$> readIO a
    set_het        a c = (\p -> c { conf_prior_het   = toProb p }) <$> readIO a

    withFd fp mode = bracket (openFd fp mode (Just 0o666) defaultFileFlags) closeFd


main :: IO ()
main = do
    (opts, files, errs) <- getOpt Permute options <$> getArgs
    unless (null errs) $ mapM_ (IO.hPutStrLn stderr) errs >> exitFailure
    conf <- foldl (>>=) (return defaultConf) opts
    case files of [f] -> main' conf f
                  _   -> IO.hPutStrLn IO.stderr "expected exactly one input file" >> exitFailure

main' :: Conf -> FilePath -> IO ()
main' Conf{..} infile = do
    conf_output $ \oiter ->
        (if infile == "-"
         then enumHandle defaultBufSize stdin
         else enumFile defaultBufSize infile) >=> run $
            joinI $ readAvroContainer $ \av_meta ->
                oiter (getRefseqs av_meta) [conf_sample] (snp_call, indel_call)
  where
    snp_call   = call conf_prior_div
    indel_call = call conf_prior_indel
    call prior lks = U.maxIndex . U.zipWith (*) lks $
                     U.replicate (U.length lks) (conf_prior_het * prior / 3)
                     U.// [(0, 1-prior)]
                     U.// [ (i, (1-conf_prior_het) * prior / 3)
                          | i <- takeWhile (< U.length lks) (scanl (+) 2 [3..]) ]



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


type CallFuncs = (U.Vector (Prob' Float) -> Int, U.Vector (Prob' Float) -> Int)
type OIter = Refs -> [S.ByteString] -> CallFuncs -> Iteratee [GenoCallBlock] IO ()
type Output a = (OIter -> IO a) -> IO a

{-
output_fasta :: FilePath -> (OIter -> IO r) -> IO r
output_fasta fn k = if fn == "-" then k (fa_out stdout)
                                 else withFile fn WriteMode $ k . fa_out
  where
    fa_out :: Handle -> Conf -> Refs -> Iteratee [Calls] IO ()
    fa_out hdl Conf{..} refs =
            by_groups p_refseq (\rs out -> do
                    let sname = sq_name $ getRef refs rs
                    out' <- lift $ enumPure1Chunk [S.concat [">", conf_sample, "--", sname]] out
                    convStream (do callz <- headStream
                                   let s1 = format_snp_call conf_prior_het callz
                                   S.append s1 <$> format_indel_call conf_prior_indel callz)
                          =$ collect_lines out') =$
            mapStreamM_ (S.hPut hdl . (flip S.snoc '\n'))
-}

-- | Formatting a SNP call.  If this was a haplopid call (four GL
-- values), we pick the most likely base and pass it on.  If it was
-- diploid, we pick the most likely dinucleotide and pass it on.
{-
format_snp_call :: Prob -> Calls -> S.ByteString
format_snp_call p cs
    | U.length gl ==  4 = S.take 1 $ S.drop (maxQualIndex gl) hapbases
    | U.length gl == 10 = S.take 1 $ S.drop (maxQualIndex $ U.zipWith (*) ps gl) dipbases
    | otherwise = error "Thou shalt not try to format_snp_call unless thou madeth a haploid or diploid call!"
  where
    gl = p_snp_pile cs
    ps = U.fromListN 10 [p,1,p,1,1,p,1,1,1,p]
    dipbases = "NAMCRSGWYKT"
    hapbases = "NACGT"
-}

-- | Formatting an Indel call.  We pick the most likely variant and
-- pass its sequence on.  Then we drop incoming calls that should be
-- deleted according to the chosen variant.  Note that this will blow up
-- unless the call was done assuming a haploid genome (which is
-- guaranteeed /in this program/)!

{-
format_indel_call :: Monad m => Prob -> Calls -> Iteratee [Calls] m S.ByteString
format_indel_call p cs
    | U.length gl0 == nv                  = go gl0
    | U.length gl0 == nv * (nv+1) `div` 2 = go homs
    | otherwise = error "Thou shalt not try to format_indel_call unless thou madeth a haploid or diploid call!"
  where
    (gl0,vars) = p_indel_pile cs
    !nv   = length vars
    !homs = U.fromListN nv [ gl0 U.! (i*(i+1) `div` 2 -1) | i <- [1..nv] ]

    go gl = I.dropWhile skip >> return (S.pack $ show $ U.toList ins)
      where
        eff_gl = U.fromList $ zipWith adjust (U.toList gl) vars
        adjust q (IndelVariant ds (V_Nuc is)) = if ds == 0 && U.null is then q else p * q

        IndelVariant del (V_Nuc ins) = ( IndelVariant 0 (V_Nuc U.empty) : vars ) !! maxQualIndex eff_gl
        skip ocs  = p_refseq ocs == p_refseq cs && p_pos ocs < p_pos cs + del

maxQualIndex :: U.Vector Prob -> Int
maxQualIndex vec = case U.ifoldl' step (0, 0, 0) vec of
    (!i, !m, !m2) -> if m / m2 > 2 then i else 0
  where
    step (!i,!m,!m2) j v = if v >= m then (j+1,v,m) else (i,m,m2)

collect_lines :: Monad m => Enumeratee S.ByteString [S.ByteString] m r
collect_lines = eneeCheckIfDone (liftI . go S.empty)
  where
    go acc k (EOF  mx) = idone (k $ Chunk [acc]) $ EOF mx
    go acc k (Chunk s) = case S.splitAt 60 (acc `S.append` s) of
                            (left, right) | S.null right -> liftI $ go left k
                                          | otherwise    -> eneeCheckIfDone (liftI . go right) . k $ Chunk [left]

-}
