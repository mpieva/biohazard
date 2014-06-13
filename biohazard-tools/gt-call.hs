{-# LANGUAGE RecordWildCards, BangPatterns, OverloadedStrings #-}
-- Command line driver for simple genotype calling.

import Bio.Adna
import Bio.Base
import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Bam.Pileup
import Bio.GenoCall
import Bio.Iteratee
import Control.Applicative
import Control.Monad
import System.Console.GetOpt
import System.Environment
import System.Exit
import System.IO

import qualified Data.ByteString.Char8          as S
import qualified Data.Iteratee                  as I
import qualified Data.Vector.Unboxed            as V

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

data Conf = Conf {
    conf_output      :: (Handle -> IO ()) -> IO (),
    conf_sample      :: S.ByteString,
    conf_ploidy      :: S.ByteString -> Int,
    conf_loverhang   :: Maybe Double,
    conf_roverhang   :: Maybe Double,
    conf_ds_deam     :: Double,
    conf_ss_deam     :: Double,
    conf_prior_het   :: Prob,
    conf_prior_indel :: Prob }

defaultConf :: Conf
defaultConf = Conf ($ stdout) "John_Doe" (const 2) Nothing Nothing
                   0.02 0.45 (qualToProb $ Q 30) (qualToProb $ Q 45)

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o" ["output"]                   (ReqArg set_output "FILE")      "Write output to FILE",
    Option "N" ["name","sample-name"]       (ReqArg set_sample "NAME")      "Set sample name to NAME",
    Option "1" ["haploid-chromosomes"]      (ReqArg set_haploid "PRF")      "Targets starting with PRF are haploid",
    Option "2" ["diploid-chromosomes"]      (ReqArg set_diploid "PRF")      "Targets starting with PRF are diploid",
    Option "l" ["overhang-length","left-overhang-length"]
                                            (ReqArg set_loverhang "LEN")    "Expected 5' overhang length is LEN",
    Option "r" ["right-overhang-length"]    (ReqArg set_roverhang "LEN")    "Expected 3' overhang length is LEN, assume single-strand prep",
    Option "d" ["deamination-rate","ds-deamination-rate","double-strand-deamination-rate"]
                                            (ReqArg set_ds_deam "FRAC")     "Deamination rate in double stranded section is FRAC",
    Option "s" ["ss-deamination-rate","single-strand-deamination-rate"]
                                            (ReqArg set_ss_deam "FRAC")     "Deamination rate in single stranded section is FRAC",
    Option "H" ["prior-heterozygous", "heterozygosity"]
                                            (ReqArg set_phet "PROB")        "Set prior for a heterozygous variant to PROB",
    -- Option "S" ["prior-snp","snp-rate","divergence"]
                                            -- (ReqArg set_pdiv "PROB")        "Set prior for an indel variant to PROB",
    Option "I" ["prior-indel","indel-rate"] (ReqArg set_pindel "PROB")      "Set prior for an indel variant to PROB",
    Option "h?" ["help","usage"]            (NoArg disp_usage)              "Display this message" ]
  where
    disp_usage _ = do pn <- getProgName
                      let blah = "Usage: " ++ pn ++ " [OPTION...] [BAM-FILE...]"
                      putStrLn $ usageInfo blah options
                      exitFailure

    set_output "-" c = return $ c { conf_output = ($ stdout) }
    set_output  fn c = return $ c { conf_output = withFile fn WriteMode }

    set_sample   nm c = return $ c { conf_sample = S.pack nm }

    set_haploid arg c = return $ c { conf_ploidy = \chr -> if S.pack arg `S.isPrefixOf` chr then 1 else conf_ploidy c chr }
    set_diploid arg c = return $ c { conf_ploidy = \chr -> if S.pack arg `S.isPrefixOf` chr then 2 else conf_ploidy c chr }

    set_loverhang a c = (\l -> c { conf_loverhang   = Just   l }) <$> readIO a
    set_roverhang a c = (\l -> c { conf_roverhang   = Just   l }) <$> readIO a
    set_ss_deam   a c = (\r -> c { conf_ss_deam     =        r }) <$> readIO a
    set_ds_deam   a c = (\r -> c { conf_ds_deam     =        r }) <$> readIO a
    set_phet      a c = (\r -> c { conf_prior_het   = toProb r }) <$> readIO a
    set_pindel    a c = (\r -> c { conf_prior_indel = toProb r }) <$> readIO a

main :: IO ()
main = do
    (opts, files, errs) <- getOpt Permute options <$> getArgs
    unless (null errs) $ mapM_ (hPutStrLn stderr) errs >> exitFailure
    Conf{..} <- foldl (>>=) (return defaultConf) opts

    let no_damage = hPutStrLn stderr "using no damage model" >> return noDamage
        ss_damage p = hPutStrLn stderr ("using single strand damage model with " ++ show p) >> return (ssDamage p)
        ds_damage p = hPutStrLn stderr ("using double strand damage model with " ++ show p) >> return (dsDamage p)

    dmg_model <- case (conf_loverhang, conf_roverhang) of
            (Nothing, Nothing) -> no_damage
            (Just ll, Nothing) -> ds_damage $ DSD conf_ss_deam conf_ds_deam (recip ll)
            (Nothing, Just lr) -> ss_damage $ SSD conf_ss_deam conf_ds_deam (recip lr) (recip lr)
            (Just ll, Just lr) -> ss_damage $ SSD conf_ss_deam conf_ds_deam (recip ll) (recip lr)

    conf_output $ \ohdl ->
        mergeInputs combineCoordinates files >=> run $ \hdr ->
            joinI $ filterStream (not . br_isUnmapped) $
            joinI $ filterStream (isValidRefseq . br_rname) $
            joinI $ by_groups same_ref (\br out -> do
                let sname = sq_name $ getRef (meta_refs hdr) $ br_rname br
                liftIO $ hPutStrLn stderr $ S.unpack sname ++ case conf_ploidy sname of
                            1 -> ": haploid call" ; 2 -> ": diploid call"
                out' <- lift $ enumPure1Chunk [S.concat [">", conf_sample, "--", sname]] out
                pileup dmg_model =$
                    mapStream (calls $! conf_ploidy sname) =$
                    convStream (headStream >>= either
                        (format_snp_call conf_prior_het)
                        (format_indel_call conf_prior_indel)) =$
                    collect_lines out') $
            mapStreamM_ (S.hPut ohdl . (flip S.snoc '\n'))

-- | This is a white lie:  We do haploid or diploid *SNP* calls, but
-- indel calls are always haploid.  Otherwise there is no good way to
-- print the result!
calls :: Int -> Pile -> Either (VarCall (GL,())) (VarCall (GL, IndelVars))
calls pl = either (Left . fmap (simple_snp_call pl)) (Right . fmap (simple_indel_call 1))

type EitherCall = Either (VarCall (GL,())) (VarCall (GL, IndelVars))

-- | Formatting a SNP call.  If this was a haplopid call (four GL
-- values), we pick the most likely base and pass it on.  If it was
-- diploid, we pick the most likely dinucleotide and pass it on.

format_snp_call :: Monad m => Prob -> VarCall (GL,()) -> Iteratee [EitherCall] m S.ByteString
format_snp_call p vc | V.length gl ==  4 = return $ S.take 1 $ S.drop (maxQualIndex gl) hapbases
                     | V.length gl == 10 = return $ S.take 1 $ S.drop (maxQualIndex $ V.zipWith (*) ps gl) dipbases
                     | otherwise = error "Thou shalt not try to format_snp_call unless thou madeth a haploid or diploid call!"
    where (gl,()) = vc_vars vc
          ps = V.fromListN 10 [p,1,p,1,1,p,1,1,1,p]
          hapbases = "NACGT"
          dipbases = "NAMCRSGWYKT"

-- | Formatting an Indel call.  We pick the most likely variant and
-- pass its sequence on.  Then we drop incoming calls that should be
-- deleted according to the chosen variant.  Note that this will blow up
-- unless the call was done assuming a haploid genome (which is
-- guaranteeed /in this program/)!
format_indel_call :: Monad m => Prob -> VarCall (GL, IndelVars) -> Iteratee [EitherCall] m S.ByteString
format_indel_call p vc | V.length gl == length vars = I.dropWhile skip >> return (S.pack $ show ins)
                       | otherwise = error "Thou shalt not have calleth format_indel_call unless thou madeth a haploid call!"
    where
        (gl,vars) = vc_vars vc
        eff_gl = V.fromList $ zipWith adjust (V.toList gl) vars
        adjust q (0,[]) = q ; adjust q _ = p * q

        (del,ins) = ( (0,[]) : vars ) !! maxQualIndex eff_gl
        skip evc  = let rs  = either vc_refseq vc_refseq evc
                        pos = either vc_pos    vc_pos    evc
                    in rs == vc_refseq vc && pos < vc_pos vc + del

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

same_ref :: BamRaw -> BamRaw -> Bool
same_ref a b = br_rname a == br_rname b

by_groups :: ( Monad m, ListLike s a, Nullable s )
          => (a -> a -> Bool) -> (a -> Enumeratee s b m r) -> Enumeratee s b m r
by_groups pr k out = do
    mhd <- peekStream
    case mhd of
        Nothing -> return out
        Just hd -> do out' <- joinI $ takeWhileE (pr hd) $ k hd out
                      by_groups pr k out'

