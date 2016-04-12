{-# LANGUAGE OverloadedStrings, RecordWildCards #-}
module Bio.Genocall.Metadata where

-- ^ Metadata necessary for a sensible genotyping workflow.

import Bio.Genocall.Adna
import Control.Applicative           hiding ( empty )
import Control.Concurrent                   ( threadDelay )
import Control.Exception
import Data.Text                            ( Text, pack )
import Data.HashMap.Strict                  ( HashMap, empty, singleton, member )
import Data.Aeson
import Data.ByteString.Char8                ( readFile, unpack )
import Data.ByteString.Lazy                 ( writeFile )
import Data.Monoid
import Data.Vector.Unboxed                  ( Vector )
import Prelude                       hiding ( writeFile, readFile )
import System.IO.Error                      ( isAlreadyExistsErrorType, ioeGetErrorType )
import System.Posix.Files.ByteString
import System.Posix.ByteString.FilePath

data Sample = Sample {
    sample_libraries   :: [Library],
    sample_avro_files  :: HashMap Text Text,                    -- ^ maps a region to the av file
    sample_bcf_files   :: HashMap Text Text,                    -- ^ maps a region to the bcf file
    sample_div_tables  :: HashMap Text (Double, Vector Int),    -- ^ maps a region to the table needed for div. estimation
    sample_divergences :: Maybe [Double]
  } deriving Show

data Library = Library {
    library_name :: Text,
    library_files :: [Text],
    library_damage :: Maybe (DamageParameters Double)
  } deriving Show

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
    toJSON (Sample ls avfs bcfs dts d) = object $ maybe id ((:) . ("divergence" .=)) d $
                                                  hashToJson "libraries" ls $
                                                  hashToJson "avro-files" avfs $
                                                  hashToJson "bcf-files" bcfs $
                                                  hashToJson "div-tables" dts []
      where
        hashToJson k vs = if null vs then id else (:) (k .= vs)

instance FromJSON Sample where
    parseJSON (String s) = pure $ Sample [Library s [s <> ".bam"] Nothing] empty empty empty Nothing
    parseJSON (Array ls) = (\ll -> Sample ll empty empty empty Nothing) <$> parseJSON (Array ls)
    parseJSON (Object o) = Sample <$> o .: "libraries"
                                  <*> (singleton "" <$> o .: "avro-file" <|> o .:? "avro-files" .!= empty)
                                  <*> (singleton "" <$> o .: "bcf-file"  <|> o .:? "bcf-files"  .!= empty)
                                  <*> o .:? "div-tables" .!= empty
                                  <*> o .:? "divergence"
    parseJSON _ = fail $ "String, Array or Object expected for Sample"

-- | Read the configuration file.  Nothing special.
readMetadata :: RawFilePath -> IO Metadata
readMetadata fn = either error return . eitherDecodeStrict =<< readFile (unpack fn)

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

split_sam_rgns :: Metadata -> [String] -> [( String, [Maybe String] )]
split_sam_rgns _meta [    ] = []
split_sam_rgns  meta (s:ss) = (s, if null rgns then [Nothing] else map Just rgns) : split_sam_rgns meta rest
    where (rgns, rest) = break (\x -> pack x `member` meta) ss


