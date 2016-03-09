{-# LANGUAGE OverloadedStrings, RecordWildCards #-}
module Bio.Genocall.Metadata where

-- ^ Metadata necessary for a sensible genotyping workflow.

import Bio.Genocall.Adna
import Control.Concurrent ( threadDelay )
import Control.Exception
import Data.Text ( Text )
import Data.HashMap.Strict ( HashMap )
import Data.Aeson
import Data.Aeson.Types
import Data.ByteString.Char8 ( readFile, unpack )
import Data.ByteString.Lazy ( writeFile )
import Data.Monoid
import Data.Vector ( fromList )
import System.IO.Error ( isAlreadyExistsErrorType, ioeGetErrorType )
import System.Posix.Files.ByteString
import System.Posix.ByteString.FilePath

import Prelude hiding ( writeFile, readFile )

data Sample = Sample {
    sample_libraries :: [Library],
    sample_divergences :: Maybe [Double] }
        deriving Show

data Library = Library {
    library_file :: Text,
    library_damage :: Maybe (DamageParameters Double) }
        deriving Show

type MetaData = HashMap Text Sample


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
    toJSON (Library file Nothing) = String file
    toJSON (Library file (Just dp)) = object [ "file" .= file, "damage" .= toJSON dp ]

instance FromJSON Library where
    parseJSON (String file) = return $ Library file Nothing
    parseJSON (Object o) = Library <$> o .: "file" <*> o .:? "damage"

instance ToJSON Sample where
    toJSON v = case v of
        Sample ls Nothing  -> toJSON ls
        Sample ls (Just d) -> object [ "libraries" .= toJSON ls , "divergence" .= toJSON d ]

instance FromJSON Sample where
    parseJSON (String s) = pure $ Sample [Library s Nothing] Nothing
    parseJSON (Array ls) = Sample <$> parseJSON (Array ls) <*> pure Nothing
    parseJSON (Object o) = Sample <$> o .: "libraries" <*> o .:? "divergence"

-- | Read the configuration file.  Nothing special.
readMetadata :: RawFilePath -> IO MetaData
readMetadata fn = either error return . eitherDecodeStrict =<< readFile (unpack fn)

-- | Update the configuration file.  First make a hard link to the
-- configuration file under a well known name (fn++"~old").  This can
-- only succeed once.  If it fails, we back off and try again later (up
-- to a very long wait).  Else, we have exclusive access.  Read the
-- file, update the data, write a new file (fn++"~new"), atomically
-- rename it and unlink the old file.
updateMetadata :: (MetaData -> MetaData) -> RawFilePath -> IO ()
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
