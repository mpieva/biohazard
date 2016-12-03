-- Genotype call a bunch of samples.  Dependency-driven, in parallel.
-- Or something.
--
-- TODO
--
-- - redeye-dar for every sample, store result
-- - redeye-pileup for every sample and every chromosome (parallel on
--   SGE), collect(!) results
-- - estimate divergence parameters
-- - run redeye-single for each sample

-- Skr1pting it is easy enough.  Can Shake do it?  I presume the
-- collection of results is somewhat less than straight forward.  Maybe
-- skr1pt first.


-- Notes:
-- - Small amounts of output land in files.  Annoying, but easy.
-- - we want -pileup and -single to run on the SGE!  (-dar and estimate
--   locally?  shouldn't be a problem)

import Bio.Bam
import Bio.Genocall.Estimators      ( estimateSingle, DivEst(..), good_regions )
import Bio.Prelude
import Data.Aeson
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

            -- final artefact: one BCF per chromosome, 'aight?
            let chromosomes = map show [1..22::Int] ++ [ "X", "Y" ]
            want [ "build/" ++ unpack (sample_name smp) ++ "." ++ chrom ++ ".bcf"
                 | chrom <- chromosomes, smp <- samples ]

            callz
            pileups
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
pileups :: Rules ()
pileups = [ "build/*.*.av", "build/*.*.divtab" ] &%> \[av,tab] -> do
                let (sm,'.':c) = splitExtension $ dropExtension $ takeFileName av
                    dmg        = "build" </> sm <.> "dmgest"
                    bams       = [ unpack libf | s <- samples, sm == unpack (sample_name s)
                                               , l <- sample_libraries s, libf <- library_files l ]
                need $ dmg : bams

                command [] "qrsh" $
                        "-now" : "no" : "-cwd" :
                        "-l" : "h_vmem=3.4G,s_vmem=3.4G,virtual_free=3.4G,s_stack=2M" :
                        "redeye-pileup" : "-o" : av : "-c" : c : "-T" : tab : "-D" : dmg : "-v" : bams

                -- command [ FileStdout tab ] "redeye-pileup" $
                --         "-o" : av : "-c" : c : "-T" : "-v" : concat libinputs


callz :: Rules ()
callz = "build/*.*.bcf" %> \bcf -> do
                let (sm,'.':c)  = splitExtension $ dropExtension $ takeFileName bcf
                    dep         = if c == "X" || c == "Y" then c else "auto"
                    divest_file = "build" </> sm <.> dep <.> "divest"
                    av_file     = "build" </> sm <.> c <.> "av"

                need [ av_file, divest_file ]

                -- this stinks.
                [dv,ht] <- either fail (return . point_est . head) . eitherDecode
                               =<< liftIO (L.readFile divest_file)

                command [] "qrsh" $
                        "-now" : "no" : "-cwd" :
                        "-l" : "h_vmem=3.4G,s_vmem=3.4G,virtual_free=3.4G,s_stack=2M" :
                        "redeye-single" : "-o" : bcf : (dropExtension bcf <.> "av")
                            : "-N" : sm : "-s"
                            : "-d" : show dv : "-D" : show ht
                            : "-i" : show (0.1*dv) : "-I" : show ht : []

-- XXX  This isn't going to work for now.
dmgests :: Rules ()
dmgests = "build/*.dmgest" %> \out -> do
                let lb = dropExtension $ takeFileName out
                    rgn_file = "build" </> lb <.> "good_regions.bam"
                need [ rgn_file ]
                command [] "qrsh" $
                        "-now" : "no" : "-cwd" :
                        "-l" : "h_vmem=3.4G,s_vmem=3.4G,virtual_free=3.4G,s_stack=2M" :
                        "redeye-dar" : {-"-T" : out :-} rgn_file : []

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
        [ "/mnt/454/HGDP/genomes_Bteam/hg19_evan.2-align/SS6004467-dedup.rg_hg19_evan.2.bam" ] ] ]


