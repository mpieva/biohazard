{-# LANGUAGE OverloadedStrings, RecordWildCards #-}
module Bio.Genocall.Metadata where

-- ^ Metadata necessary for a sensible genotyping workflow.

import Bio.Genocall.Adna
import Control.Applicative
import Control.Concurrent                   ( threadDelay )
import Control.Exception
import Data.Text                            ( Text )
import Data.HashMap.Strict                  ( HashMap, traverseWithKey )
import Data.Aeson
import Data.ByteString.Char8                ( readFile, unpack )
import Data.ByteString.Lazy                 ( writeFile )
import Data.Monoid
import Prelude                       hiding ( writeFile, readFile )
import System.IO.Error                      ( isAlreadyExistsErrorType, ioeGetErrorType )
import System.Posix.Files.ByteString
import System.Posix.ByteString.FilePath

data Sample = Sample {
    sample_libraries :: [Library],
    sample_avro_file :: Text,
    sample_bcf_file :: Text,
    sample_divergences :: Maybe [Double] }
        deriving Show

data Library = Library {
    library_name :: Text,
    library_files :: [Text],
    library_damage :: Maybe (DamageParameters Double) }
        deriving Show

newtype Metadata' = Metadata' Metadata
type Metadata = HashMap Text Sample


instance ToJSON float => ToJSON (DamageParameters float) where
    toJSON DP{..} = object [ "ss-sigma"  .= ssd_sigma
                           , "ss-delta"  .= ssd_delta
                           , "ss-lambda" .= ssd_lambda
                           , "ss-kappa"  .= ssd_kappa
                           , "ds-sigma"  .= dsd_sigma
                           , "ds-delta"  .= dsd_delta
                           , "ds-lambda" .= dsd_lambda ]

instance FromJSON float => FromJSON (DamageParameters float) where
    parseJSON = withObject "damage parameters" $ \o ->
                    DP <$> o .: "ss-sigma"
                       <*> o .: "ss-delta"
                       <*> o .: "ss-lambda"
                       <*> o .: "ss-kappa"
                       <*> o .: "ds-sigma"
                       <*> o .: "ds-delta"
                       <*> o .: "ds-lambda"

instance ToJSON Library where
    toJSON (Library name files dp) = object ( maybe id ((:) . ("damage" .=)) dp
                                            $ [ "name" .= name, "files" .= files ] )

instance FromJSON Library where
    parseJSON (String name) = return $ Library name [name <> ".bam"] Nothing
    parseJSON (Object o) = Library <$> o .: "name"
                                   <*> (maybe id (:) <$> o .:? "file"
                                                     <*> o .:? "files" .!= [])
                                   <*> o .:? "damage"
    parseJSON _ = fail "String or Object expected for library"

instance ToJSON Sample where
    toJSON (Sample ls avf bcf d) = object $ maybe id ((:) . ("divergence" .=)) d $
                                          [ "libraries" .= ls, "avro-file" .= avf, "bcf-file"  .= bcf ]

instance FromJSON Metadata' where
    parseJSON = withObject "metadata must be an object" $ fmap Metadata' . traverseWithKey p_sample
      where
        p_sample nm (String s) = pure $ Sample [Library s [s <> ".bam"] Nothing] (nm <> ".av") (nm <> ".bcf") Nothing
        p_sample nm (Array ls) = (\ll -> Sample ll (nm <> ".av") (nm <> ".bcf") Nothing) <$> parseJSON (Array ls)
        p_sample nm (Object o) = Sample <$> o .: "libraries"
                                        <*> o .:? "avro-file" .!= (nm <> ".av")
                                        <*> o .:? "bcf-file" .!= (nm <> ".bcf")
                                        <*> o .:? "divergence"
        p_sample nm _ = fail $ "String, Array or Object expected for Sample " ++ show nm

-- | Read the configuration file.  Nothing special.
readMetadata :: RawFilePath -> IO Metadata
readMetadata fn = either error (\(Metadata' m) -> return m) . eitherDecodeStrict =<< readFile (unpack fn)

-- | Update the configuration file.  First make a hard link to the
-- configuration file under a well known name (fn++"~old").  This can
-- only succeed once.  If it fails, we back off and try again later (up
-- to a very long wait).  Else, we have exclusive access.  Read the
-- file, update the data, write a new file (fn++"~new"), atomically
-- rename it and unlink the old file.
updateMetadata :: (Metadata -> Metadata) -> RawFilePath -> IO ()
updateMetadata f fp = go (36::Int)     -- retry every 5 secs for 3 minutes
  where
    go n = catchJust
                (\e -> if isAlreadyExistsErrorType (ioeGetErrorType e) && n > 0 then Just () else Nothing)
                (do createLink fp (fp <> "~old")
                    mdata <- readMetadata fp
                    writeFile (unpack $ fp <> "~new") . encode . toJSON $ f mdata
                    rename (fp <> "~new") fp
                `finally` removeLink (fp <> "~old"))
                (\_ -> threadDelay 5000000 >> go (n-1))
