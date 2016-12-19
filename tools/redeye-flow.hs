-- Genotype call a bunch of samples.  Dependency-driven, in parallel.
-- Or something.
--
-- TODO
--
-- - redeye-dar for every sample, store result
-- - redeye-single for every sample and every chromosome
-- - estimate divergence for each sample

-- Skr1pting it is easy enough.  Can Shake do it?  I presume the
-- collection of results is somewhat less than straight forward.  Maybe
-- skr1pt first.

-- Notes:
-- - Small amounts of output land in files.  Annoying, but easy.
-- - We run -dar and -single on the SGE.  Other parts run locally.

import Bio.Bam
import Bio.Genocall.Estimators      ( estimateSingle, good_regions )
import Bio.Prelude
import Data.Aeson.Encode.Pretty
import Data.Binary                  ( decodeOrFail )
import Development.Shake
import Development.Shake.FilePath
import System.Directory
import System.IO

import qualified Data.ByteString.Lazy   as L
import qualified Data.Sequence          as Z
import qualified Data.Vector.Generic    as V

data Sample = Sample {
    sample_name      :: Text,
    sample_libraries :: [ Library ]
  } deriving Show

data Library = Library {
    library_name :: Text,
    library_files :: [ Text ]
  } deriving Show



main :: IO ()
main = shakeArgs shakeOptions { shakeFiles = "_shake" } $ do

            -- final artefacts: one BCF per chromosome,
            let chromosomes = map show [1..22::Int] ++ [ "X", "Y" ]
            want [ "build/" ++ unpack (sample_name smp) ++ "." ++ chrom ++ ".bcf"
                 | chrom <- chromosomes, smp <- samples ]

            -- and div/het estimates
            want [ "build/" ++ unpack (sample_name smp) ++ "." ++ part ++ ".divest"
                 | part <- ["auto","X","Y"], smp <- samples ]

            callz
            divests
            dmgests
            rgn_files


divests :: Rules ()
divests = do
            "build/*.auto.divest" %> \out -> do
                let stem = dropExtension $ dropExtension out
                lReadFiles' [ stem ++ "." ++ show c ++ ".divtab" | c <- [1..22::Int] ] >>=
                    either fail_decode (\tabs -> liftIO $ do
                            (de1,de2) <- estimateSingle $ mconcat [ t | (_,_,t) <- tabs ]
                            L.writeFile out $ encodePretty [ de1, de2 ])
                        . mapM decodeOrFail

            "build/*.X.divest" %> \out -> do
                lReadFile' (out -<.> "divtab") >>=
                    either fail_decode (\(_,_,tab) -> liftIO $ do
                            (de1,de2) <- estimateSingle tab
                            L.writeFile out $ encodePretty [ de1, de2 ])
                        . decodeOrFail

            "build/*.Y.divest" %> \out -> do
                lReadFile' (out -<.> "divtab") >>=
                    either fail_decode (\(_,_,tab) -> liftIO $ do
                            (de1,de2) <- estimateSingle tab
                            L.writeFile out $ encodePretty [ de1, de2 ])
                        . decodeOrFail
  where
    fail_decode (rest,off,msg) = error $
        msg ++ " at " ++ shows off " near " ++ show (L.take 16 rest)

    lReadFile'   x = need [x] >> liftIO (L.readFile x)
    lReadFiles' xs = need  xs >> liftIO (mapM L.readFile xs)


-- one pileup per chrmosome * sample; input is the
-- bam files and one dmgest per sample
callz :: Rules ()
callz = [ "build/*.*.bcf", "build/*.*.divtab" ] &%> \[bcf,tab] -> do
                let (sm,'.':c) = splitExtension $ dropExtension $ takeFileName bcf
                    dmg        = "build" </> sm <.> "dmgest"
                    bams       = [ unpack libf | s <- samples, sm == unpack (sample_name s)
                                               , l <- sample_libraries s, libf <- library_files l ]
                need $ dmg : bams

                command [] "qrsh" $
                        "-now" : "no" : "-cwd" :
                        "-l" : "h_vmem=3.4G,s_vmem=3.4G,virtual_free=3.4G,s_stack=2M" :
                        "redeye-single" : "-o" : bcf : "-c" : c : "-T" : tab : "-D" : dmg
                                        : "-N" : sm : "-v" : bams

                -- command [ FileStdout tab ] "redeye-pileup" $
                --         "-o" : av : "-c" : c : "-T" : "-v" : concat libinputs


dmgests :: Rules ()
dmgests = "build/*.dmgest" %> \out -> do
                let lb = dropExtension $ takeFileName out
                    rgn_file = "build" </> lb <.> "good_regions.bam"
                need [ rgn_file ]
                command [] "qrsh" $
                        "-now" : "no" : "-cwd" :
                        "-l" : "h_vmem=3.4G,s_vmem=3.4G,virtual_free=3.4G,s_stack=2M" :
                        "redeye-dar" : "-o" : out : rgn_file : []

rgn_files :: Rules ()
rgn_files = "build/*.good_regions.bam" %> \out -> do
                let sm = dropExtension $ dropExtension $ takeFileName out
                    lfs = [ unpack f | s <- samples, unpack (sample_name s) == sm
                                     , l <- sample_libraries s
                                     , f <- library_files l ]
                need lfs

                liftIO $ subsetbams out lfs (takeLen 1000000 good_regions) 35
  where
    takeLen !n (x@(_,_,l):xs) | n > 0 = x : takeLen (n-l) xs
    takeLen  _             _          = []


-- | Reads regions from many bam files, writes one.
-- XXX  It might make sense to serialize not BAM, but the result of
-- piling up.
subsetbams :: FilePath -> [FilePath] -> [( Bytes, Int, Int )] -> Int -> IO ()
subsetbams ofp (ifp:ifps) rgns0 minlen = do
    withFile (ofp ++ "~") WriteMode                             $ \hdl ->
        go ifp ifps >=> run                                         $ \hdr ->
        filterStream ((>= minlen) . V.length . b_seq . unpackBam)  =$
        writeBamHandle hdl hdr
    renameFile (ofp ++ "~") ofp
  where
    enum1 :: (MonadIO m, MonadMask m) => FilePath -> Enumerator' BamMeta [BamRaw] m a
    enum1 fp k = do idx <- liftIO $ readBamIndex fp
                    enumFileRandom defaultBufSize fp >=> run >=> run $
                        decodeAnyBam $ \hdr ->
                            let rgns = sort [ Region (Refseq $ fromIntegral ri) p (p+l)
                                            | (ch, p, l) <- rgns0
                                            , let Just ri = Z.findIndexL ((==) ch . sq_name) (meta_refs hdr) ]
                            in eneeBamRegions idx rgns (k hdr)

    go :: (MonadIO m, MonadMask m) => FilePath -> [FilePath] -> Enumerator' BamMeta [BamRaw] m b
    go fp [       ] = enum1 fp
    go fp (fp1:fps) = mergeEnums' (go fp1 fps) (enum1 fp) combineCoordinates

subsetbams ofp [] rgns0 minlen =
    error $ "Wait, what? " ++ show (ofp, rgns0, minlen)

samples :: [Sample]
samples =
    [ let lib nm = Library nm [ nm <> ".bam" ]
      in  Sample "HC" (map lib [ "A9368", "A9369", "A9401", "A9402", "A9403", "A9404", "B8747", "R5473" ])
    , let lane i = Library (fromString $ show (i::Int)) [ fromString (printf path i) ]
          path   = "/mnt/ngs_data/140411_SN7001204_0257_AC2MW7ACXX_PEdi_SP/Ibis/BWA/proc1/s_%d_sequence_ancient_hg19_evan.bam"
      in Sample "Vanity" (map lane [3..8])
    , Sample "SS6004467" [ Library "SS6004467"
        [ "/mnt/454/HGDP/genomes_Bteam/hg19_evan.2-align/SS6004467-dedup.rg_hg19_evan.2.bam" ] ]

    , Sample "Goyet"
        [ Library "A9122" [ "/mnt/expressions/mateja/Goyet/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9122_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9229" [ "/mnt/expressions/mateja/Goyet/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9229_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9349" [ "/mnt/expressions/mateja/Goyet/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9349_final_sorted.uniq.L35MQ0.bam" ] ]

    , Sample "Vindija_G1"
        [ Library "A9121" [ "/mnt/expressions/mateja/Vindija_G1/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9121_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9228" [ "/mnt/expressions/mateja/Vindija_G1/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9228_final_sorted.uniq.L35MQ0.bam" ] ]
        -- , Library "A9248" [ "/mnt/expressions/mateja/Vindija_G1/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9248_final_sorted.uniq.L35MQ0.bam" ] ]

    , Sample "Les_Cottes"
        [ Library "A9230" [ "/mnt/expressions/mateja/Les_Cottes/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9230_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9290" [ "/mnt/expressions/mateja/Les_Cottes/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9290_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9291" [ "/mnt/expressions/mateja/Les_Cottes/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9291_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9309" [ "/mnt/expressions/mateja/Les_Cottes/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9309_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9350" [ "/mnt/expressions/mateja/Les_Cottes/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9350_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9393" [ "/mnt/expressions/mateja/Les_Cottes/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9393_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9394" [ "/mnt/expressions/mateja/Les_Cottes/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9394_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9395" [ "/mnt/expressions/mateja/Les_Cottes/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9395_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9420" [ "/mnt/expressions/mateja/Les_Cottes/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9420_final_sorted.uniq.L35MQ0.bam" ] ]

    , Sample "Mezmaiskaya2"
        [ Library "A9180" [ "/mnt/expressions/mateja/Mezmaiskaya2/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9180_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9288" [ "/mnt/expressions/mateja/Mezmaiskaya2/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9288_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9289" [ "/mnt/expressions/mateja/Mezmaiskaya2/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9289_final_sorted.uniq.L35MQ0.bam" ] ]

    , Sample "Spy"
        [ Library "A9416" [ "/mnt/expressions/mateja/Spy/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9416_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9417" [ "/mnt/expressions/mateja/Spy/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9417_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9418" [ "/mnt/expressions/mateja/Spy/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9418_final_sorted.uniq.L35MQ0.bam" ]
        , Library "A9419" [ "/mnt/expressions/mateja/Spy/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/A9419_final_sorted.uniq.L35MQ0.bam" ]
        , Library "R5556" [ "/mnt/expressions/mateja/Spy/FinalBAMs/Per_library/AnalyzeBAM_L35MQ0_per_library/R5556_final_sorted.uniq.L35MQ0.bam" ] ] ]
