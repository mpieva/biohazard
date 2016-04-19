{-# LANGUAGE OverloadedStrings, RecordWildCards, FlexibleContexts, BangPatterns #-}
module Bio.Genocall.Metadata where

-- ^ Metadata necessary for a sensible genotyping workflow.

import Bio.Genocall.Adna                    ( DamageParameters(..) )
import Control.Applicative           hiding ( empty )
import Control.Concurrent                   ( threadDelay )
import Control.Exception                    ( bracket, onException, handleJust )
import Control.Monad                        ( forM_ )
import Data.Text                            ( Text, pack )
import Data.HashMap.Strict                  ( HashMap )
import Data.Aeson
import Data.ByteString.Char8                ( readFile )
import Data.ByteString.Lazy                 ( toChunks )
import Data.ByteString.Unsafe               ( unsafeUseAsCStringLen )
import Data.Monoid
import Data.Vector.Unboxed                  ( Vector )
import Foreign.Ptr                          ( castPtr )
import GHC.IO.Exception                     ( IOErrorType(..) )
import Prelude                       hiding ( writeFile, readFile )
import System.IO.Error                      ( isAlreadyExistsErrorType, ioeGetErrorType )
import System.Posix.Files                   ( rename, removeLink )
import System.Posix.IO

import qualified Data.HashMap.Strict as M

data Sample = Sample {
    sample_libraries   :: [Library],
    sample_avro_files  :: HashMap Text Text,                    -- ^ maps a region to the av file
    sample_bcf_files   :: HashMap Text Text,                    -- ^ maps a region to the bcf file
    sample_div_tables  :: HashMap Text (Double, Vector Int),    -- ^ maps a region to the table needed for div. estimation
    sample_divergences :: HashMap Text DivEst
  } deriving Show

data Library = Library {
    library_name :: Text,
    library_files :: [Text],
    library_damage :: Maybe (DamageParameters Double)
  } deriving Show

-- | Divergence estimate.  Lists contain three or four floats, these are
-- divergence, heterozygosity at W sites, heterozygosity at S sites, and
-- optionally gappiness in this order.
data DivEst = DivEst {
    point_est :: [Double],
    conf_region :: [( [Double], [Double] )]
  } deriving Show


type Metadata = HashMap Text Sample

instance ToJSON DivEst where
    toJSON DivEst{..} = object $ [ "estimate" .= point_est
                                 , "confidence-region" .= conf_region ]

instance FromJSON DivEst where
    parseJSON (Object o) = DivEst <$> o .: "estimate" <*> o .:? "confidence-region" .!= []
    parseJSON (Array a) = flip DivEst [] <$> parseJSON (Array a)
    parseJSON _ = fail $ "divergence estimate should be an array or an object"


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
    toJSON (Sample ls avfs bcfs dts ds) = object $ hashToJson "divergences" ds   $
                                                   listToJson "libraries"   ls   $
                                                   hashToJson "avro-files"  avfs $
                                                   hashToJson "bcf-files"   bcfs $
                                                   hashToJson "div-tables"  dts  []
      where
        hashToJson k vs = if M.null vs then id else (:) (k .= vs)
        listToJson k vs = if   null vs then id else (:) (k .= vs)

instance FromJSON Sample where
    parseJSON (String s) = pure $ Sample [Library s [s <> ".bam"] Nothing] M.empty M.empty M.empty M.empty
    parseJSON (Array ls) = (\ll -> Sample ll M.empty M.empty M.empty M.empty) <$> parseJSON (Array ls)
    parseJSON (Object o) = Sample <$> o .: "libraries"
                                  <*> (M.singleton "" <$> o .: "avro-file" <|> o .:? "avro-files" .!= M.empty)
                                  <*> (M.singleton "" <$> o .: "bcf-file"  <|> o .:? "bcf-files"  .!= M.empty)
                                  <*> o .:? "div-tables" .!= M.empty
                                  <*> (M.singleton "" <$> o .: "divergence" <|> o.:? "divergences" .!= M.empty)
    parseJSON _ = fail $ "String, Array or Object expected for Sample"


-- | Read the configuration file.  Retries, because NFS tends to result
-- in 'ResourceVanished' if the file is replaced while we try to read it.
readMetadata :: FilePath -> IO Metadata
readMetadata fn = either error return . eitherDecodeStrict =<< go (15::Int)
  where
    go !n = handleJust     -- retry every sec for 15 seconds
                (\e -> case ioeGetErrorType e of ResourceVanished | n > 0 -> Just () ; _ -> Nothing)
                (\_ -> threadDelay 1000000 >> go (n-1))
                (readFile fn)

-- | Update the configuration file.  Open a new file (fn++"~new") in
-- exclusive mode.  Then read the old file, write the update to the new
-- file, rename it atomically, then close it.  Use of O_EXCL should
-- ensure that nobody interferes.  This is atomic even on NFS, provided
-- NFS and kernel are new enough.  For older NFS, I cannot be bothered.
--
-- (The first idea was to base this on the supposed fact that link(2) is
-- atomic and fails if the new filename exists.  This approach does seem
-- to contain a race condition, though.)
updateMetadata :: (Metadata -> Metadata) -> FilePath -> IO ()
updateMetadata f fp = go (36::Int)     -- retry every 5 secs for 3 minutes
  where
    fpn = fp <> "~new"

    go !n = handleJust
                (\e -> if isAlreadyExistsErrorType (ioeGetErrorType e) && n > 0 then Just () else Nothing)
                (\_ -> threadDelay 5000000 >> go (n-1)) $ do
                bracket
                    (openFd fpn WriteOnly (Just 0o666) defaultFileFlags{ exclusive = True })
                    (closeFd) $ \fd ->
                        (do mdata <- readMetadata fp
                            forM_ (toChunks . encode . toJSON $ f mdata) $ \ch ->
                                unsafeUseAsCStringLen ch $ \(p,l) ->
                                    fdWriteBuf fd (castPtr p) (fromIntegral l)
                            rename fpn fp)
                        `onException` removeLink fpn


split_sam_rgns :: Metadata -> [String] -> [( String, [Maybe String] )]
split_sam_rgns _meta [    ] = []
split_sam_rgns  meta (s:ss) = (s, if null rgns then [Nothing] else map Just rgns) : split_sam_rgns meta rest
    where (rgns, rest) = break (\x -> pack x `M.member` meta) ss
