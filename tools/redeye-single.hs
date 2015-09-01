{-# LANGUAGE BangPatterns, RecordWildCards, OverloadedStrings #-}
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
-- circumstances, else VCF or BCF would be canonical.  Or maybe VCF
-- restricted to variant sites.  Or VCF restricted to sites not known to
-- be reference.  People will come up with filters for sure...

import Bio.Base
import Bio.Bam
import Bio.Bam.Pileup
import Bio.Genocall.AvroFile
import Bio.Iteratee
import Control.Monad
import Data.Avro
import Data.Foldable ( toList )
import System.Console.GetOpt
import System.Environment
import System.Exit
import System.IO

import qualified Data.ByteString.Char8          as S
import qualified Data.Vector.Unboxed            as U

data Conf = Conf {
    conf_output      :: Maybe (Output ()),
    conf_sample      :: S.ByteString,
    conf_ploidy      :: S.ByteString -> Int,
    conf_report      :: String -> IO (),
    conf_prior_div   :: Prob Double,
    conf_prior_het   :: Prob Double,
    conf_prior_indel :: Prob Double }

defaultConf :: Conf
defaultConf = Conf Nothing "John_Doe" (const 2) (\_ -> return ())
                   (qualToProb $ Q 30) (qualToProb $ Q 20) (qualToProb $ Q 45)

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

    be_verbose       c = return $ c { conf_report = hPutStrLn stderr }
    set_sample     a c = return $ c { conf_sample = S.pack a }
    set_bcf_out  "-" c = return $ c { conf_output = undefined }
    set_bcf_out   fn c = return $ c { conf_output = undefined }

    set_hap a c = return $ c { conf_ploidy = \chr -> if S.pack a `S.isPrefixOf` chr then 1 else conf_ploidy c chr }
    set_dip a c = return $ c { conf_ploidy = \chr -> if S.pack a `S.isPrefixOf` chr then 2 else conf_ploidy c chr }

    set_div        a c = (\p -> c { conf_prior_div   = toProb p }) <$> readIO a
    set_indel      a c = (\p -> c { conf_prior_indel = toProb p }) <$> readIO a
    set_het        a c = (\p -> c { conf_prior_het   = toProb p }) <$> readIO a


main :: IO ()
main = do
    (opts, files, errs) <- getOpt Permute options <$> getArgs
    unless (null errs) $ mapM_ (hPutStrLn stderr) errs >> exitFailure
    conf <- foldl (>>=) (return defaultConf) opts
    case files of [f] -> main' conf f
                  _   -> hPutStrLn stderr "expected exactly one input file" >> exitFailure

main' :: Conf -> FilePath -> IO ()
main' Conf{..} infile = do
    -- Pieces of old code were no good at all...  so, reboot.
    -- We read an AVRO container, we apply the priors, we format the
    -- result as minimal BCF (because VCF is an abomination).
    --
    -- So, generate a header (VCF or BCF) first.  We get the /names/ of
    -- the reference sequences from the schema and the /lengths/ from
    -- the biohazard.refseq_length entry.
    --
    -- Then app priors to each record, call genotype and format.
    --
    -- Need to write BCF2.  Which is custom binary goop inside a BGZF
    -- container.  :(
    enumFile defaultBufSize infile >=> run $
        joinI $ readAvroContainer $ \av_meta -> do
            -- needs /names/ and /lengths/ of reference sequences
            -- names come from the enum definition, need to find the schema
            -- lengths are stuck in biohazard.reference_lengths
            liftIO $ S.hPut stdout $
                vcf_header (getRefseqs av_meta)
            a <- peekStream
            return $ const () (a :: Maybe GenoCallBlock)



vcf_header :: Refs -> S.ByteString
vcf_header refs = S.unlines . map (S.append "##") $
    [ "fileformat=VCFv4.2"
    , "INFO=<ID=MQ,Number=1,Type=Integer,Description=\"RMS mapping quality\">"
    , "INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Number of MAPQ==0 reads covering this record\">"
    , "FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    , "FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"read depth\">"
    , "FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"genotype likelihoods in deciban\">"
    , "FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"conditional genotype quality in deciban\">" ] ++
    [ S.concat [ "contig=<ID=", sq_name s, ",length=", S.pack (show (sq_length s)), ">" ] | s <- toList refs ]


type OIter = Conf -> Refs -> Iteratee [Calls] IO ()
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
format_snp_call :: Prob Double -> Calls -> S.ByteString
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
format_indel_call :: Monad m => Prob Double -> Calls -> Iteratee [Calls] m S.ByteString
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
-}

maxQualIndex :: U.Vector (Prob Double) -> Int
maxQualIndex vec = case U.ifoldl' step (0, 0, 0) vec of
    (!i, !m, !m2) -> if m / m2 > 2 then i else 0
  where
    step (!i,!m,!m2) j v = if v >= m then (j+1,v,m) else (i,m,m2)

{-
collect_lines :: Monad m => Enumeratee S.ByteString [S.ByteString] m r
collect_lines = eneeCheckIfDone (liftI . go S.empty)
  where
    go acc k (EOF  mx) = idone (k $ Chunk [acc]) $ EOF mx
    go acc k (Chunk s) = case S.splitAt 60 (acc `S.append` s) of
                            (left, right) | S.null right -> liftI $ go left k
                                          | otherwise    -> eneeCheckIfDone (liftI . go right) . k $ Chunk [left]

-}
