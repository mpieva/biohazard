{-# LANGUAGE CPP, DeriveGeneric #-}
-- | Metadata necessary for a sensible genotyping workflow.
module Bio.Genocall.Metadata where

import Bio.Adna                             ( DamageParameters(..), NewDamageParameters(..) )
import Bio.Prelude                   hiding ( writeFile, readFile )
import Bio.Genocall.Estimators
import Bio.Util.Pretty
import Data.Aeson
import Data.Binary
import Data.Binary.Get                      ( runGetOrFail )
import Data.Binary.Put                      ( runPut )
import Data.ByteString.Lazy                 ( readFile )
import Data.Text.Encoding                   ( decodeUtf8', encodeUtf8 )
import System.Posix.Files                   ( rename, removeLink )
import System.Posix.IO

import qualified Data.HashMap.Strict as M
import qualified Data.Vector.Unboxed as U

data Sample = Sample {
    sample_libraries   :: [Library],
    sample_avro_files  :: M.HashMap Text Text,                    -- ^ maps a region to the av file
    sample_bcf_files   :: M.HashMap Text Text,                    -- ^ maps a region to the bcf file
    sample_div_tables  :: M.HashMap Text DivTable,                -- ^ maps a region to the table needed for div. estimation
    sample_divergences :: M.HashMap Text DivEst
  } deriving Show

data Library = Library {
    library_name :: Text,
    library_files :: [Text],
    library_damage :: GenDamageParameters U.Vector Double
  } deriving Show

type Metadata = M.HashMap Text Sample

instance ToJSON DivEst where
    toJSON DivEst{..} = object $ [ "estimate" .= point_est
                                 , "confidence-region" .= conf_region ]

instance Binary DivEst where
    put DivEst{..} = put point_est >> put conf_region
    get = DivEst <$> get <*> get

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

instance (ToJSON (vec float), ToJSON float) => ToJSON (NewDamageParameters vec float) where
    toJSON NDP{..} = object [ "gc-frac" .= dp_gc_frac
                            , "mu"      .= dp_mu
                            , "nu"      .= dp_nu
                            , "alpha5"  .= dp_alpha5
                            , "beta5"   .= dp_beta5
                            , "alpha"   .= dp_alpha
                            , "beta"    .= dp_beta
                            , "alpha3"  .= dp_alpha3
                            , "beta3"   .= dp_beta3 ]

instance Binary float => Binary (DamageParameters float) where
    put DP{..} = put ssd_sigma >> put ssd_delta >> put ssd_lambda >> put ssd_kappa >>
                 put dsd_sigma >> put dsd_delta >> put dsd_lambda
    get = DP <$> get <*> get <*> get <*> get <*> get <*> get <*> get

instance (Binary float, U.Unbox float) => Binary (NewDamageParameters U.Vector float) where
    put NDP{..} = put dp_gc_frac >> put dp_mu >> put dp_nu >>
                  putVector dp_alpha5 >> putVector dp_beta5  >>
                  put dp_alpha >> put dp_beta >>
                  putVector dp_alpha3 >> putVector dp_beta3
    get = NDP <$> get <*> get <*> get <*> getVector <*> getVector
              <*> get <*> get <*> getVector <*> getVector

instance (Binary float, U.Unbox float) => Binary (GenDamageParameters U.Vector float) where
    put  UnknownDamage  = put (0::Word8)
    put (OldDamage  dp) = put (1::Word8) >> put dp
    put (NewDamage  dp) = put (2::Word8) >> put dp

    get = get >>= \x -> case x :: Word8 of 0 -> pure UnknownDamage
                                           1 -> OldDamage <$> get
                                           2 -> NewDamage <$> get


instance FromJSON float => FromJSON (DamageParameters float) where
    parseJSON = withObject "damage parameters" $ \o ->
                    DP <$> o .: "ss-sigma"
                       <*> o .: "ss-delta"
                       <*> o .: "ss-lambda"
                       <*> o .: "ss-kappa"
                       <*> o .: "ds-sigma"
                       <*> o .: "ds-delta"
                       <*> o .: "ds-lambda"

instance (FromJSON (vec float), FromJSON float) => FromJSON (NewDamageParameters vec float) where
    parseJSON = withObject "new damage parameters" $ \o ->
                    NDP <$> o .: "gc-frac"
                        <*> o .: "mu"
                        <*> o .: "nu"
                        <*> o .: "alpha5"
                        <*> o .: "beta5"
                        <*> o .: "alpha"
                        <*> o .: "beta"
                        <*> o .: "alpha3"
                        <*> o .: "beta3"

instance ToJSON Library where
    toJSON (Library name files dp) = object $ ( case dp of
                                                  UnknownDamage -> id
                                                  OldDamage dp' -> (:) ("damage" .= dp')
                                                  NewDamage dp' -> (:) ("damage" .= dp')
                                              ) [ "name" .= name, "files" .= files ]

instance Binary Library where
    put Library{..} = put library_name >> put library_files >> put library_damage
    get = Library <$> get <*> get <*> get

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

instance Binary Sample where
    put Sample{..} = put sample_libraries >> putObject sample_avro_files >>
                     putObject sample_bcf_files >> putObject sample_div_tables >>
                     putObject sample_divergences
    get = Sample <$> get <*> getObject <*> getObject <*> getObject <*> getObject

instance FromJSON Sample where
    parseJSON (String s) = pure $ Sample [Library s [s <> ".bam"] UnknownDamage] M.empty M.empty M.empty M.empty
    parseJSON (Array ls) = (\ll -> Sample ll M.empty M.empty M.empty M.empty) <$> parseJSON (Array ls)
    parseJSON (Object o) = Sample <$> o .: "libraries"
                                  <*> (M.singleton "" <$> o .: "avro-file" <|> o .:? "avro-files" .!= M.empty)
                                  <*> (M.singleton "" <$> o .: "bcf-file"  <|> o .:? "bcf-files"  .!= M.empty)
                                  <*> o .:? "div-tables" .!= M.empty
                                  <*> (M.singleton "" <$> o .: "divergence" <|> o.:? "divergences" .!= M.empty)
    parseJSON _ = fail $ "String, Array or Object expected for Sample"

instance FromJSON DivTable where
    parseJSON x = parseJSON x >>= \[a,b] -> DivTable <$> parseJSON a <*> parseJSON b

instance Binary DivTable where
    put (DivTable a b) = put a >> putVector b
    get = DivTable <$> get <*> getVector

instance ToJSON DivTable where
    toJSON (DivTable a b) = toJSON [toJSON a, toJSON b]

#if !MIN_VERSION_text(1,2,1)
instance Binary Text where
    put = put . encodeUtf8
    get = do
      bs <- get
      case decodeUtf8' bs of
          Left exn -> fail $ show exn
          Right  a -> return a
#endif

putObject :: Binary value => M.HashMap Text value -> Put
putObject m = put (M.size m) >> M.foldrWithKey (\k v a -> put k >> put v >> a) (return ()) m

getObject :: Binary value => Get (M.HashMap Text value)
getObject = get >>= \l -> get_map M.empty (l::Int)
    where
        get_map !acc 0 = return acc
        get_map !acc n = get >>= \k -> get >>= \v -> get_map (M.insert k v acc) (n-1)

-- Hm.  I'm putting the vector in reverse order, because the
-- accumulation when reading reverses it.
putVector :: (U.Unbox a, Binary a) => U.Vector a -> Put
putVector v = put (U.length v) >> U.mapM_ put (U.reverse v)

getVector :: (U.Unbox a, Binary a) => Get (U.Vector a)
getVector = get >>= \l -> U.fromListN l <$> get_list [] l
    where
        get_list acc 0 = return acc
        get_list acc n = get >>= \ !x -> get_list (x:acc) (n-1)


-- | Read the configuration file.  Retries, because NFS tends to result
-- in 'ResourceVanished' if the file is replaced while we try to read it.
readJsonMetadata :: FilePath -> IO Metadata
readJsonMetadata fn = either error return . eitherDecode =<< go (15::Int)
  where
    go !n = handleJust     -- retry every sec for 15 seconds
                (\e -> case ioeGetErrorType e of ResourceVanished | n > 0 -> Just () ; _ -> Nothing)
                (\_ -> threadDelay 1000000 >> go (n-1))
                (readFile fn)

-- | Read the configuration file.  Retries, because NFS tends to result
-- in 'ResourceVanished' if the file is replaced while we try to read it.
readMetadata :: FilePath -> IO Metadata
readMetadata fn = either (error . show) (\(_,_,r) -> return r) . runGetOrFail getObject =<< go (15::Int)
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
updateMetadata f fp = go (360::Int)     -- retry every 5 secs for 30 minutes
  where
    fpn = fp <> "~new"

    go !n = handleJust
                (\e -> if isAlreadyExistsErrorType (ioeGetErrorType e) && n > 0 then Just () else Nothing)
                (\_ -> threadDelay 5000000 >> go (n-1)) $ do
                bracket
                    (openFd fpn WriteOnly (Just 0o666) defaultFileFlags{ exclusive = True })
                    (closeFd) $ \fd ->
                        (do mdata <- readMetadata fp
                            fdPutLazy fd . runPut . putObject $ f mdata
                            rename fpn fp)
                        `onException` removeLink fpn

writeMetadata :: FilePath -> Metadata -> IO ()
writeMetadata fp mdata = bracket
                    (openFd fp WriteOnly (Just 0o666) defaultFileFlags{ exclusive = True })
                    (closeFd) $ \fd ->
                        fdPutLazy fd . runPut . putObject $ mdata

split_sam_rgns :: Metadata -> [String] -> [( String, [Maybe String] )]
split_sam_rgns _meta [    ] = []
split_sam_rgns  meta (s:ss) = (s, if null rgns then [Nothing] else map Just rgns) : split_sam_rgns meta rest
    where (rgns, rest) = break (\x -> fromString x `M.member` meta) ss
