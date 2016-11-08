-- Genotype call a bunch of samples.  Dependency-driven, in parallel.
-- Or something.
--
-- TODO
--
-- - redeye-dar for every library
--   * store result
-- - redeye-pileup for every sample and every chromosome (parallel on SGE)
--   * collect(!) results
-- - estimate divergence parameters
-- - run redeye-single for each sample

-- Skr1pting it is easy enough.  Can Shake do it?  I presume the
-- collection of results is somewhat less than straight forward.  Maybe
-- skr1pt first.


-- Notes:
-- - Small amounts of output land in files.  Annoying, but easy.
-- - we want -pileup and -single to run on the SGE!  (-dar and estimate
--   locally?  shouldn't be a problem)

import Bio.Genocall.Estimators      ( estimateSingle, DivEst(..) )
import Bio.Prelude
import Bio.Util.Pretty              ( pshow, pparse )
import Data.Binary                  ( decodeOrFail )

import Development.Shake
import Development.Shake.FilePath

import qualified Data.ByteString.Lazy   as L
import qualified Data.Text.IO           as S
import qualified Data.Text.Lazy         as T
import qualified Data.Text.Lazy.IO      as T

data Sample = Sample {
    sample_name      :: Text,
    sample_libraries :: [ Library ]
  } deriving Show

data Library = Library {
    library_name :: Text,
    library_files :: [ Text ]
  } deriving Show


{-
instance FromJSON Library where
    parseJSON (String name) = return $ Library name [name <> ".bam"] UnknownDamage
    parseJSON (Object o) = Library <$> o .: "name"
                                   <*> (maybe id (:) <$> o .:? "file"
                                                     <*> o .:? "files" .!= [])
                                   <*> (OldDamage <$> o .: "damage" <|>
                                        NewDamage <$> o .: "damage" <|>
                                        pure UnknownDamage)
    parseJSON _ = fail "String or Object expected for library"

instance ToJSON Sample where
    toJSON (Sample ls avfs bcfs dts ds) = object $ hashToJson "divergences" ds   $
                                                   listToJson "libraries"   ls   $
                                                   hashToJson "avro-files"  avfs $
                                                   hashToJson "bcf-files"   bcfs $
                                                   hashToJson "div-tables"  dts  []
      where
        hashToJson k vs = if M.null vs then id else (:) (k .= vs)
        listToJson k vs = if   null vs then id else (:) (k .= vs)

instance FromJSON Sample where
    parseJSON (String s) = pure $ Sample [Library s [s <> ".bam"] UnknownDamage] M.empty M.empty M.empty M.empty
    parseJSON (Array ls) = (\ll -> Sample ll M.empty M.empty M.empty M.empty) <$> parseJSON (Array ls)
    parseJSON (Object o) = Sample <$> o .: "libraries"
                                  <*> (M.singleton "" <$> o .: "avro-file" <|> o .:? "avro-files" .!= M.empty)
                                  <*> (M.singleton "" <$> o .: "bcf-file"  <|> o .:? "bcf-files"  .!= M.empty)
                                  <*> o .:? "div-tables" .!= M.empty
                                  <*> (M.singleton "" <$> o .: "divergence" <|> o.:? "divergences" .!= M.empty)
    parseJSON _ = fail $ "String, Array or Object expected for Sample"
-}

main :: IO ()
main = shakeArgs shakeOptions { shakeFiles = "_shake" } $ do

            mapM_ buildSample samples

            "build/*.auto.divest" %> \out -> do
                let stem = dropExtension $ dropExtension out
                lReadFiles' [ stem ++ "." ++ show c ++ ".divtab" | c <- [1..22::Int] ] >>=
                    either fail_decode (\tabs -> liftIO $ do
                            (de1,de2) <- estimateSingle $ mconcat [ t | (_,_,t) <- tabs ]
                            T.writeFile out $ T.unlines [ pshow de1, pshow de2 ])
                        . mapM decodeOrFail

            "build/*.X.divest" %> \out -> do
                lReadFile' (out -<.> "divtab") >>=
                    either fail_decode (\(_,_,tab) -> liftIO $ do
                            (de1,de2) <- estimateSingle tab
                            T.writeFile out $ T.unlines [ pshow de1, pshow de2 ])
                        . decodeOrFail

            "build/*.Y.divest" %> \out -> do
                lReadFile' (out -<.> "divtab") >>=
                    either fail_decode (\(_,_,tab) -> liftIO $ do
                            (de1,de2) <- estimateSingle tab
                            T.writeFile out $ T.unlines [ pshow de1, pshow de2 ])
                        . decodeOrFail
  where
    fail_decode (rest,off,msg) = error $
        msg ++ " at " ++ shows off " near " ++ show (L.take 16 rest)

    lReadFile'   x = need [x] >> liftIO (L.readFile x)
    lReadFiles' xs = need  xs >> liftIO (mapM L.readFile xs)


buildSample :: Sample -> Rules ()
buildSample smp = do
    -- final artefact: one BCF per chromosome, 'aight?
    let chromosomes = map show [1..22::Int] ++ [ "X", "Y" ]
    want [ "build/" ++ unpack (sample_name smp) ++ "." ++ chrom ++ ".bcf" | chrom <- chromosomes ]

    -- want [ "build/" ++ unpack (sample_name smp) ++ "." ++ subset ++ ".divest"
         -- | subset <- [ "auto", "X", "Y" ] ]

    [ "build/" ++ unpack (sample_name smp) ++ ".*.av",
      "build/" ++ unpack (sample_name smp) ++ ".*.divtab" ] &%> \[av,tab] -> do
        need [ "build/" ++ unpack (library_name lib) ++ ".dmgest" | lib <- sample_libraries smp ]
        let '.':c = takeExtension $ dropExtension av

        libinputs <- sequence [ do de <- readFile' $ "build/" ++ unpack (library_name lib) ++ ".dmgest"
                                   return $ "-D" : flatten de : map unpack (library_files lib)
                              | lib <- sample_libraries smp ]

        command [] "qrsh" $
                "-now" : "no" : "-cwd" :
                "-l" : "h_vmem=3.4G,s_vmem=3.4G,virtual_free=3.4G,s_stack=2M" :
                "redeye-pileup" : "-o" : av : "-c" : c : "-T" : tab : "-v" : concat libinputs

        -- command [ FileStdout tab ] "redeye-pileup" $
        --         "-o" : av : "-c" : c : "-T" : "-v" : concat libinputs

    "build/" ++ unpack (sample_name smp) ++ ".*.bcf" %> \bcf -> do
        let '.':c = takeExtension $ dropExtension bcf
            dep = if c == "X" || c == "Y" then c else "auto"
            divest_file = dropExtension (dropExtension bcf) ++ "." ++ dep ++ ".divest"

        need [ dropExtension bcf <.> "av", divest_file ]

        [dv,ht] <- either fail (return . point_est) . pparse =<< liftIO (S.readFile divest_file)

        command [] "qrsh" $
                "-now" : "no" : "-cwd" :
                "-l" : "h_vmem=3.4G,s_vmem=3.4G,virtual_free=3.4G,s_stack=2M" :
                "redeye-single" : "-o" : bcf : (dropExtension bcf <.> "av")
                    : "-N" : unpack (sample_name smp) : "-s"
                    : "-d" : show dv : "-D" : show ht
                    : "-i" : show (0.1*dv) : "-I" : show ht : []


    mapM_ buildLib $ sample_libraries smp
  where
    flatten = map $ \c -> if c == '\n' then ' ' else c

buildLib :: Library -> Rules ()
buildLib lib = do
    "build/" ++ unpack (library_name lib) ++ ".dmgest" %> \out -> do
        let lfs = map unpack $ library_files lib
        need lfs
        command [ FileStdout out ] "redeye-dar" ("-m" : "35" : lfs )

samples :: [Sample]
samples = [ Sample "HC" (map lib [ "A9368", "A9369", "A9401", "A9402", "A9403", "A9404", "B8747", "R5473" ]) ]
  where
    lib nm = Library nm [ nm <> ".bam" ]

