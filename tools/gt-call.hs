{-# LANGUAGE RecordWildCards, BangPatterns, OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell, FlexibleContexts #-}
-- Command line driver for simple genotype calling.

import Bio.Base
import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Bam.Pileup
import Bio.Genocall
import Bio.Genocall.Adna
import Bio.Genocall.AvroFile
import Bio.Iteratee
import Bio.Util                                 ( float2mini )
import Control.Applicative
import Control.DeepSeq
import Control.Monad
import Data.Avro
import Data.Function
import System.Console.GetOpt
import System.Environment
import System.Exit
import System.IO

-- import qualified Data.ByteString                as B
import qualified Data.ByteString.Char8          as S
import qualified Data.Iteratee                  as I
-- import qualified Data.Text                      as T
import qualified Data.Text.Encoding             as T
import qualified Data.Vector.Unboxed            as V

-- import Debug.Trace

-- Ultimately, we might produce a VCF file looking somewhat like this:
--
-- ##FORMAT=<ID=A,Number=2,Type=Integer,Description="Number of A bases on forward and reverse strand">
-- ##FORMAT=<ID=C,Number=2,Type=Integer,Description="Number of C bases on forward and reverse strand">
-- ##FORMAT=<ID=G,Number=2,Type=Integer,Description="Number of G bases on forward and reverse strand">
-- ##FORMAT=<ID=T,Number=2,Type=Integer,Description="Number of T bases on forward and reverse strand">
--      (we should count bases on both strands for this)
--
-- ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (only filtered reads used for calling)">
-- ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
-- ##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
--      (basic statistics. we keep these)
--
-- ##FORMAT=<ID=IR,Number=1,Type=Integer,Description="Number of reads with InDel starting at this position">
-- ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
-- ##INFO=<ID=Dels,Number=1,Type=Float,Description="Fraction of Reads Containing Spanning Deletions">
--      (this is bullshit)
--
-- ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
-- ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
-- ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
--      (these are straight forward to compute?)
--
-- ##INFO=<ID=AF1000g,Number=1,Type=Float,Description="Global alternative allele frequency (AF)...">
-- ##INFO=<ID=AMR_AF,Number=1,Type=Float,Description="Alternative allele frequency (AF) for samples from AMR based on 1000G">
-- ##INFO=<ID=ASN_AF,Number=1,Type=Float,Description="Alternative allele frequency (AF) for samples from ASN based on 1000G">
-- ##INFO=<ID=AFR_AF,Number=1,Type=Float,Description="Alternative allele frequency (AF) for samples from AFR based on 1000G">
-- ##INFO=<ID=EUR_AF,Number=1,Type=Float,Description="Alternative allele frequency (AF) for samples from EUR based on 1000G">
-- ##INFO=<ID=1000gALT,Number=1,Type=String,Description="Alternative allele referred to by 1000G">
-- ##INFO=<ID=TS,Number=1,Type=String,Description="Sequences in Ensembl v64 EPO Compara 6 primate block">
-- ##INFO=<ID=TSseq,Number=1,Type=String,Description="Primary species bases (in order of TS field) in the EPO Compara 6 primate block">
-- ##INFO=<ID=CAnc,Number=1,Type=String,Description="Ref-Chimp/Human ancestor base at this position">
-- ##INFO=<ID=GAnc,Number=1,Type=String,Description="Ref-Gorilla ancestor base at this position">
-- ##INFO=<ID=OAnc,Number=1,Type=String,Description="Ref-Orang ancestor base at this position">
-- ##INFO=<ID=mSC,Number=1,Type=Float,Description="PhastCons Mammalian conservation score (excluding human)">
-- ##INFO=<ID=pSC,Number=1,Type=Float,Description="PhastCons Primate conservation score (excluding human)">
-- ##INFO=<ID=GRP,Number=1,Type=Float,Description="GERP conservation score">
-- ##INFO=<ID=bSC,Number=1,Type=Float,Description="B score">
-- ##INFO=<ID=Map20,Number=1,Type=Float,Description="Mapability score of Duke University (determined from 20bp reads)">
-- ##INFO=<ID=RM,Number=0,Type=Flag,Description="Position is repeat masked in the reference sequence of the EPO 6 primate block">
-- ##INFO=<ID=SysErr,Number=0,Type=Flag,Description="Position was identified as systematic error in the 1000 genome trios">
-- ##INFO=<ID=SysErrHCB,Number=0,Type=Flag,Description="Position was identified as systematic error based on shared SNPs...">
-- ##INFO=<ID=UR,Number=0,Type=Flag,Description="Position is in a copy number control region identified by the Eichler lab">
--      (this is external, will not be generated)
--
-- ##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position is in a CpG context based on the Ref/Ancestor">
-- ##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods...">
--      (this is computable, isn't it?!)
--
-- ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
--      (this is from VarScan 2, a program that uses fixed cutoffs.  It
--      is not clear that this has any use at all.)
--
-- ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
-- ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
-- ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
-- ##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
-- ##INFO=<ID=DP,Number=1,Type=Integer,Description="Filtered Depth">
-- ##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
-- ##INFO=<ID=HRun,Number=1,Type=Integer,Description="Largest Contiguous Homopolymer Run of Variant Allele In Either Direction">
-- ##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
-- ##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
-- ##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
-- ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
--      (WTF?)

-- parameters used for the Unified Genotyper:
--      downsample_to_coverage=250
--      heterozygosity=0.001
--      pcr_error_rate=1.0E-4
--      indel_heterozygosity=1.25E-4


-- auxilliary files (from Martin's option parser):
--
--      ancestor_path       EMF     /mnt/expressions/martin/sequence_db/epo/epo_6_primate_v64/split/
--      G1000               VCF     /mnt/expressions/martin/sequence_db/snps/20110521_G1000_release/phase1_intergrated_calls.20101123.snps_indels_svs.sites.vcf.gz
--      bscores             TSV1i   /mnt/454/Altaiensis/users/martin/HighCoverage/additional_information/bscores/liftover/human.tsv.gz
--      mammalscores        TSV2f   /mnt/454/Altaiensis/users/martin/HighCoverage/additional_information/mammal_conservation/liftover/human.tsv.gz
--      primatescores       TSV2f   /mnt/454/Altaiensis/users/martin/HighCoverage/additional_information/primate_conservation/liftover/human.tsv.gz
--      gerpscores          TSV2f   /mnt/454/Altaiensis/users/fernando/sequencedb/GERP/liftover/human.tsv.gz
--      mapability          TSV2i   /mnt/454/Altaiensis/users/martin/HighCoverage/additional_information/mapability/liftover/human.tsv.gz
--      uregions            TSV1    /mnt/454/Altaiensis/users/martin/HighCoverage/additional_information/EL_control_regions/liftover/human.tsv.gz
--      syserrors           TSV1    /mnt/454/Altaiensis/users/martin/HighCoverage/additional_information/sys_errors/liftover/human.tsv.gz
--      syserrorsHCB        TSV1    /mnt/454/Altaiensis/users/fernando/sequencedb/SysErrHCB/human.tsv.gz

--  TSV1:  chr start end score
--  TSV2:  chr pos score

-- About damage parameters:  We effectively have three different models
-- (SS, DS, no damage) and it may not be possible to choose one a
-- priori.  To manage this cleanly, we should have one universal model,
-- but the three we have are not generalizations of each other.
--
-- So we treat the choice of model as another parameter.  We feed
-- parameters for all three in, together with probabilities for each.
-- Said probabilities are derived from the likelihoods obtained when
-- fitting the parameters individually.  Genotype calling then involves
-- calling once under each model and summing them (effectively
-- marginalizing on the choice of model).

data Conf = Conf {
    conf_output      :: Maybe Output,
    conf_sample      :: S.ByteString,
    conf_ploidy      :: S.ByteString -> Int,
    conf_damage      :: Maybe (DamageParameters Double),
    conf_loverhang   :: Maybe Double,
    conf_roverhang   :: Maybe Double,
    conf_ds_deam     :: Double,
    conf_ss_deam     :: Double,
    conf_theta       :: Maybe Double,
    conf_report      :: String -> IO (),
    conf_prior_het   :: Prob,
    conf_prior_indel :: Prob }

defaultConf :: Conf
defaultConf = Conf Nothing "John_Doe" (const 2) Nothing Nothing Nothing
                   0.02 0.45 Nothing (\_ -> return ())
                   (qualToProb $ Q 30) (qualToProb $ Q 45)

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o" ["output", "avro-output"]    (ReqArg set_avro_out "FILE")    "Write AVRO output to FILE",
    Option [ ] ["fasta-output"]             (ReqArg set_fa_output "FILE")   "Write FA output to FILE",
    Option "N" ["name","sample-name"]       (ReqArg set_sample "NAME")      "Set sample name to NAME",
    Option "1" ["haploid-chromosomes"]      (ReqArg set_haploid "PRF")      "Targets starting with PRF are haploid",
    Option "2" ["diploid-chromosomes"]      (ReqArg set_diploid "PRF")      "Targets starting with PRF are diploid",
    Option "D" ["damage"]                   (ReqArg set_damage "PARMS")     "Set universal damage parameters",
    Option "l" ["overhang-param","left-overhang-param"]
                                            (ReqArg set_loverhang "PROB")   "Parameter for 5' overhang length is PROB",
    Option "r" ["right-overhang-param"]     (ReqArg set_roverhang "PROB")   "Parameter for 3' overhang length is PROB, assume single-strand prep",
    Option "d" ["deamination-rate","ds-deamination-rate","double-strand-deamination-rate"]
                                            (ReqArg set_ds_deam "FRAC")     "Deamination rate in double stranded section is FRAC",
    Option "s" ["ss-deamination-rate","single-strand-deamination-rate"]
                                            (ReqArg set_ss_deam "FRAC")     "Deamination rate in single stranded section is FRAC",
    Option "t" ["theta","dependency-coefficient"]
                                            (ReqArg set_theta   "FRAC")     "Set dependency coefficient to FRAC (\"N\" to turn off)",
    Option "H" ["prior-heterozygous", "heterozygosity"]
                                            (ReqArg set_phet "PROB")        "Set prior for a heterozygous variant to PROB",
    -- Removed this, because it needs access to a reference.
    -- But maybe we can derive this from a suitable BAM file?
    -- Or move it to another tool?
    -- Option "S" ["prior-snp","snp-rate","divergence"]
                                            -- (ReqArg set_pdiv "PROB")        "Set prior for an indel variant to PROB",
    Option "I" ["prior-indel","indel-rate"] (ReqArg set_pindel "PROB")      "Set prior for an indel variant to PROB",
    Option "v" ["verbose"]                  (NoArg be_verbose)              "Print more diagnostics",
    Option "h?" ["help","usage"]            (NoArg disp_usage)              "Display this message" ]
  where
    disp_usage _ = do pn <- getProgName
                      let blah = "Usage: " ++ pn ++ " [OPTION...] [BAM-FILE...]"
                      putStrLn $ usageInfo blah options
                      exitFailure

    be_verbose c = return $ c { conf_report = hPutStrLn stderr }

    set_fa_output fn = add_output $ output_fasta fn
    set_avro_out  fn = add_output $ output_avro  fn

    add_output ofn cf =
        return $ cf { conf_output = Just $ \k ->
            ofn $ \oit1 -> maybe (k oit1) ($ \oit2 -> k (\c r -> () <$ I.zip (oit1 c r) (oit2 c r))) (conf_output cf) }

    set_sample   nm c = return $ c { conf_sample = S.pack nm }

    set_haploid arg c = return $ c { conf_ploidy = \chr -> if S.pack arg `S.isPrefixOf` chr then 1 else conf_ploidy c chr }
    set_diploid arg c = return $ c { conf_ploidy = \chr -> if S.pack arg `S.isPrefixOf` chr then 2 else conf_ploidy c chr }

    set_theta "N" c = return $ c { conf_theta       =  Nothing }
    set_theta     a c = (\t -> c { conf_theta       = Just   t }) <$> readIO a
    set_loverhang a c = (\l -> c { conf_loverhang   = Just   l }) <$> readIO a
    set_roverhang a c = (\l -> c { conf_roverhang   = Just   l }) <$> readIO a
    set_ss_deam   a c = (\r -> c { conf_ss_deam     =        r }) <$> readIO a
    set_ds_deam   a c = (\r -> c { conf_ds_deam     =        r }) <$> readIO a
    set_phet      a c = (\r -> c { conf_prior_het   = toProb r }) <$> readIO a
    set_pindel    a c = (\r -> c { conf_prior_indel = toProb r }) <$> readIO a
    set_damage    a c = (\u -> c { conf_damage      = Just   u }) <$> readIO a

main :: IO ()
main = do
    (opts, files, errs) <- getOpt Permute options <$> getArgs
    unless (null errs) $ mapM_ (hPutStrLn stderr) errs >> exitFailure
    conf@Conf{..} <- foldl (>>=) (return defaultConf) opts

    let no_damage   = conf_report "using no damage model" >> return noDamage
        ss_damage p = conf_report ("using single strand damage model with " ++ show p) >> return (univDamage p)
        ds_damage p = conf_report ("using double strand damage model with " ++ show p) >> return (univDamage p)
        u_damage  p = conf_report ("using universal damage parameters " ++ show p) >> return (univDamage p)

    dmg_model <- case (conf_damage, conf_loverhang, conf_roverhang) of
            (Just u,        _, _) -> u_damage u
            (_, Nothing, Nothing) -> no_damage
            (_, Just pl, Nothing) -> ds_damage $ DP 0 0 0 0 conf_ss_deam conf_ds_deam pl
            (_, Nothing, Just pr) -> ss_damage $ DP conf_ss_deam conf_ds_deam pr pr 0 0 0
            (_, Just pl, Just pr) -> ss_damage $ DP conf_ss_deam conf_ds_deam pl pr 0 0 0

    maybe (output_fasta "-") id conf_output $ \oiter ->
        mergeInputs combineCoordinates files >=> run $ \hdr ->
            filterStream (not . br_isUnmapped) =$
            filterStream (isValidRefseq . br_rname) =$
            progressPos "GT call at " conf_report (meta_refs hdr) =$
            by_groups ((==) `on` br_rname) (\br out -> do
                let sname = sq_name $ getRef (meta_refs hdr) $ br_rname br
                    pl = conf_ploidy sname
                liftIO $ conf_report $ S.unpack sname ++ ["",": haploid call",": diploid call"] !! pl
                pileup dmg_model =$ mapStream (calls conf_theta pl) out) =$
            oiter conf (meta_refs hdr)


type OIter = Conf -> Refs -> Iteratee [Calls] IO ()
type Output = (OIter -> IO ()) -> IO ()

output_fasta :: FilePath -> (OIter -> IO r) -> IO r
output_fasta fn k = if fn == "-" then k (fa_out stdout)
                                 else withFile fn WriteMode $ k . fa_out
  where
    fa_out :: Handle -> Conf -> Refs -> Iteratee [Calls] IO ()
    fa_out hdl Conf{..} refs =
            by_groups ((==) `on` p_refseq) (\cs out -> do
                    let sname = sq_name $ getRef refs $ p_refseq cs
                    out' <- lift $ enumPure1Chunk [S.concat [">", conf_sample, "--", sname]] out
                    convStream (do callz <- headStream
                                   let s1 = format_snp_call conf_prior_het callz
                                   S.append s1 <$> format_indel_call conf_prior_indel callz)
                          =$ collect_lines out') =$
            mapStreamM_ (S.hPut hdl . (flip S.snoc '\n'))


-- | We do calls of any ploidy, but the FastA output code will fail if
-- the ploidy isn't 1 or 2.  For indel calls, the FastA output will also
-- cheat and pretend it was a haploid call.
--
-- XXX  For the time being, forward and reverse piles get concatenated.
-- For the naive call, this doesn't matter.  For the MAQ call, it feels
-- more correct to treat them separately and multiply (add?) the results.

calls :: Maybe Double -> Int -> Pile -> Calls
calls Nothing pl pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !s = simple_snp_call pl $ uncurry (++) $ p_snp_pile pile
    !i = force $ simple_indel_call pl $ p_indel_pile pile

calls (Just theta) pl pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !s = maq_snp_call pl theta $ uncurry (++) $ p_snp_pile pile -- XXX
    !i = force $ simple_indel_call pl $ p_indel_pile pile

-- | Formatting a SNP call.  If this was a haplopid call (four GL
-- values), we pick the most likely base and pass it on.  If it was
-- diploid, we pick the most likely dinucleotide and pass it on.

format_snp_call :: Prob -> Calls -> S.ByteString
format_snp_call p cs
    | V.length gl ==  4 = S.take 1 $ S.drop (maxQualIndex gl) hapbases
    | V.length gl == 10 = S.take 1 $ S.drop (maxQualIndex $ V.zipWith (*) ps gl) dipbases
    | otherwise = error "Thou shalt not try to format_snp_call unless thou madeth a haploid or diploid call!"
  where
    gl = p_snp_pile cs
    ps = V.fromListN 10 [p,1,p,1,1,p,1,1,1,p]
    dipbases = "NAMCRSGWYKT"
    hapbases = "NACGT"

-- | Formatting an Indel call.  We pick the most likely variant and
-- pass its sequence on.  Then we drop incoming calls that should be
-- deleted according to the chosen variant.  Note that this will blow up
-- unless the call was done assuming a haploid genome (which is
-- guaranteeed /in this program/)!

format_indel_call :: Monad m => Prob -> Calls -> Iteratee [Calls] m S.ByteString
format_indel_call p cs
    | V.length gl0 == nv                  = go gl0
    | V.length gl0 == nv * (nv+1) `div` 2 = go homs
    | otherwise = error "Thou shalt not try to format_indel_call unless thou madeth a haploid or diploid call!"
  where
    (gl0,vars) = p_indel_pile cs
    !nv   = length vars
    !homs = V.fromListN nv [ gl0 V.! (i*(i+1) `div` 2 -1) | i <- [1..nv] ]

    go gl = I.dropWhile skip >> return (S.pack $ show $ V.toList ins)
      where
        eff_gl = V.fromList $ zipWith adjust (V.toList gl) vars
        adjust q (IndelVariant ds (V_Nuc is)) = if ds == 0 && V.null is then q else p * q

        IndelVariant del (V_Nuc ins) = ( IndelVariant 0 (V_Nuc V.empty) : vars ) !! maxQualIndex eff_gl
        skip ocs  = p_refseq ocs == p_refseq cs && p_pos ocs < p_pos cs + del

maxQualIndex :: V.Vector Prob -> Int
maxQualIndex vec = case V.ifoldl' step (0, 0, 0) vec of
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

by_groups :: ( Monad m, ListLike s a, Nullable s )
          => (a -> a -> Bool) -> (a -> Enumeratee s b m r) -> Enumeratee s b m r
by_groups pr k out = do
    mhd <- peekStream
    case mhd of
        Nothing -> return out
        Just hd -> takeWhileE (pr hd) =$ k hd out >>= by_groups pr k


output_avro :: FilePath -> (OIter -> IO r) -> IO r
output_avro fn k = if fn == "-" then k (av_out stdout)
                                else withFile fn WriteMode $ k . av_out
  where
    av_out :: Handle -> Conf -> Refs -> Iteratee [Calls] IO ()
    av_out hdl _cfg refs = compileBlocks refs =$
                           writeAvroContainer ContainerOpts{..} =$
                           mapChunksM_ (S.hPut hdl)

    objects_per_block = 16
    filetype_label = "Genotype Likelihoods V0.1"


-- Serialize the results from genotype calling in a sensible way.  We
-- write an Avro file, but we add another blocking layer on top so we
-- don't need to endlessly repeat coordinates.

compileBlocks :: Monad m => Refs -> Enumeratee [Calls] [GenoCallBlock] m a
compileBlocks refs = convStream $ do
        c1 <- headStream
        tailBlock (p_refseq c1) (p_pos c1) (p_pos c1) (16*1024 :: Int) [pack c1]
  where
    tailBlock !rs !p0 !po !n acc = do
        mc <- peekStream
        case mc of
            Just c1 | rs == p_refseq c1 && po+1 == p_pos c1 && n > 0 -> do
                    _ <- headStream
                    tailBlock rs p0 (po+1) (n-1) $ pack c1 : acc

            _ -> return [ GenoCallBlock
                    { reference_name = T.decodeLatin1 $ sq_name $ getRef refs rs
                    , start_position = p0
                    , called_sites   = reverse acc } ]

    pack c1 = GenoCallSite{..}
      where
        snp_stats         = p_snp_stat c1
        indel_stats       = p_indel_stat c1
        snp_likelihoods   = compact_likelihoods $ p_snp_pile c1
        indel_likelihoods = compact_likelihoods $ fst $ p_indel_pile c1
        indel_variants    = snd $ p_indel_pile c1

-- | Storing likelihoods:  we take the natural logarithm (GL values are
-- already in a log scale) and convert to minifloat 0.4.4
-- representation.  Range and precision should be plenty.
compact_likelihoods :: V.Vector Prob -> V.Vector Word8
compact_likelihoods = V.map (float2mini . negate . unPr)

