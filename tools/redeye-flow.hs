-- Genotype call a bunch of samples.  Dependency-driven, in parallel.
--
-- - Small amounts of output land in files.  Annoying, but easy.
-- - We optionally run -dar and -single on the SGE.  Other parts run locally.

import Bio.Bam
import Bio.Genocall.Estimators      ( estimateSingle, good_regions )
import Bio.Prelude
import Data.Aeson
import Data.Aeson.Encode.Pretty
import Data.Aeson.Types
import Data.Binary                  ( decodeOrFail )
import Development.Shake
import Development.Shake.FilePath
import System.Directory
import System.Console.GetOpt
import System.IO

import qualified Data.ByteString        as B
import qualified Data.ByteString.Lazy   as L
import qualified Data.HashMap.Strict    as H
import qualified Data.Sequence          as Z
import qualified Data.Vector.Generic    as V

data Sample = Sample {
    sample_name      :: Text,
    sample_libraries :: [ Text ]
  } deriving Show

parseSamples :: Value -> Parser [Sample]
parseSamples = withObject "samples" $ \o ->
    sequence [ Sample k <$> parseStrings v | (k,v) <- H.toList o ]
  where
    parseStrings v = withText "file name" (return . (:[])) v
                 <|> withArray "file names"
                     (mapM (withText "file name" return) . V.toList) v


data Flags = GridEngine | Samples FilePath deriving Eq

flagOptions :: [ OptDescr (Either String Flags) ]
flagOptions = [ Option "G" ["grid-engine"] (NoArg $ Right GridEngine) "Run on grid engine (uses qrsh)."
              , Option "S" ["samples"] (ReqArg (Right . Samples) "FILE") "Read samples from FILE" ]


main :: IO ()
main = shakeArgsWith shakeOptions flagOptions $ \flags targets -> return $ Just $ do
            samples <- liftIO $ concat <$> sequence
                        [ do str <- B.readFile fp
                             case eitherDecodeStrict' str of
                                Left err -> error err
                                Right val -> case parseEither parseSamples val of
                                    Left err -> error err
                                    Right smps -> return smps
                        | Samples fp <- flags ]

            if null targets then do
                -- final artefacts: one BCF per chromosome,
                let chromosomes = map show [1..22::Int] ++ [ "X", "Y" ]
                want [ "build/" ++ unpack (sample_name smp) ++ "." ++ chrom ++ ".bcf"
                     | chrom <- chromosomes, smp <- samples ]

                -- and div/het estimates
                want [ "build/" ++ unpack (sample_name smp) ++ "." ++ part ++ ".divest"
                     | part <- ["auto","X","Y"], smp <- samples ]
              else
                want targets

            callz samples flags
            divests
            dmgests flags
            rgn_files samples


divests :: Rules ()
divests = do
            "build/*.auto.divest" %> \out -> do
                let stem = dropExtension $ dropExtension out
                raw <- lReadFiles' [ stem ++ "." ++ show c ++ ".divtab" | c <- [1..22::Int] ]
                putLoud $ "estimate for " ++ out
                either fail_decode (\tabs -> liftIO $ do
                            (de1,de2) <- estimateSingle $ mconcat [ t | (_,_,t) <- tabs ]
                            L.writeFile out $ encodePretty [ de1, de2 ])
                        $ mapM decodeOrFail raw

            "build/*.X.divest" %> \out -> do
                raw <- lReadFile' (out -<.> "divtab")
                putLoud $ "estimate for " ++ out
                either fail_decode (\(_,_,tab) -> liftIO $ do
                            (de1,de2) <- estimateSingle tab
                            L.writeFile out $ encodePretty [ de1, de2 ])
                        $ decodeOrFail raw

            "build/*.Y.divest" %> \out -> do
                raw <- lReadFile' (out -<.> "divtab")
                putLoud $ "estimate for " ++ out
                either fail_decode (\(_,_,tab) -> liftIO $ do
                            (de1,de2) <- estimateSingle tab
                            L.writeFile out $ encodePretty [ de1, de2 ])
                        $ decodeOrFail raw
  where
    fail_decode (rest,off,msg) = error $
        msg ++ " at " ++ shows off " near " ++ show (L.take 16 rest)

    lReadFile'   x = need [x] >> liftIO (L.readFile x)
    lReadFiles' xs = need  xs >> liftIO (mapM L.readFile xs)


-- one pileup per chrmosome * sample; input is the
-- bam files and one dmgest per sample
callz :: [Sample] -> [Flags] -> Rules ()
callz samples flags = [ "build/*.*.bcf", "build/*.*.divtab" ] &%> \[bcf,tab] -> do
                let (sm,'.':c) = splitExtension $ dropExtension $ takeFileName bcf
                    dmg        = "build" </> sm <.> "dmgest"
                    bams       = [ unpack libf | s <- samples, sm == unpack (sample_name s)
                                               , libf <- sample_libraries s ]
                need $ dmg : bams

                if GridEngine `elem` flags
                    then
                        unsafeExtraThread $
                            command [] "qrsh" $
                                "-now" : "no" : "-cwd" : "-N" : (sm ++ "-" ++ c) :
                                "-l" : "h_vmem=3.4G,s_vmem=3.4G,virtual_free=3.4G,s_stack=2M" :
                                "redeye-single" : "-o" : bcf : "-c" : c : "-T" : tab : "-D" : dmg
                                                : "-N" : sm : "-v" : bams
                    else
                        command [] "redeye-single" $
                            "-o" : bcf : "-c" : c : "-T" : tab : "-D" : dmg
                                 : "-N" : sm : "-v" : bams


dmgests :: [Flags] -> Rules ()
dmgests flags = "build/*.dmgest" %> \out -> do
                let sm = dropExtension $ takeFileName out
                    rgn_file = "build" </> sm <.> "good_regions.bam"
                need [ rgn_file ]

                if GridEngine `elem` flags
                    then
                        unsafeExtraThread $
                            command [] "qrsh" $
                                "-now" : "no" : "-cwd" : "-N" : (sm ++ "-dar") :
                                "-l" : "h_vmem=3.4G,s_vmem=3.4G,virtual_free=3.4G,s_stack=2M" :
                                "redeye-dar" : "-o" : out : rgn_file : []
                    else
                        command [] "redeye-dar" $ "-o" : out : rgn_file : []

rgn_files :: [Sample] -> Rules ()
rgn_files samples = do
    mem <- newResource "heavy IO" 1
    "build/*.good_regions.bam" %> \out -> do
                let sm = dropExtension $ dropExtension $ takeFileName out
                    lfs = [ unpack f | s <- samples, unpack (sample_name s) == sm
                                     , f <- sample_libraries s ]
                need lfs
                withResource mem 1 $ do
                    putLoud $ "subsetting " ++ show lfs
                    liftIO $ subsetbams out lfs (takeLen 5000000 good_regions) 35
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


