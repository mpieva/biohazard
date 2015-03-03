{-# LANGUAGE OverloadedStrings, FlexibleInstances, TemplateHaskell #-}
{-# LANGUAGE RecordWildCards, BangPatterns, FlexibleContexts #-}
module Data.Avro where

import Bio.Iteratee
import Control.Applicative
import Control.Monad
import Control.Monad.ST ( runST, ST )
import Data.Aeson hiding ((.=))
import Data.Array.MArray
import Data.Array.ST ( STUArray )
import Data.Array.Unsafe ( castSTUArray )
import Data.Binary.Get
import Data.Bits
import Data.Binary.Builder
import Data.Foldable ( foldMap )
import Data.Int ( Int64 )
import Data.Maybe
import Data.Monoid
import Data.Scientific
import Data.Text.Encoding
import Data.Word ( Word32, Word64 )
import Foreign.Storable ( Storable, sizeOf )
import Language.Haskell.TH
import System.Random

import qualified Data.ByteString as B
import qualified Data.ByteString.Lazy as BL
import qualified Data.HashMap.Strict as H
import qualified Data.ListLike as LL
import qualified Data.Text as T
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U

-- ^ Support for Avro.
-- Current status is that we can generate schemas for certain Haskell
-- values, serialize to binary and JSON representations, and write
-- Container files using the null codec.  The C implementation likes
-- some, but not all of these containers; it's unclear if that's the
-- fault of the C implementation, though.
--
-- Meanwhile, serialization works for nested sums-of-products, as long as the
-- product uses record syntax and the top level is a plain record.
-- The obvious primitives are supported.

(.=) :: ToJSON a => String -> a -> (T.Text, Value)
k .= v = (T.pack k, toJSON v)

string :: String -> Value
string = String . T.pack

-- | This is the class of types we can embed into the Avro
-- infrastructure.  Right now, we can derive a schema, encode to
-- the Avro binary format, and encode to the Avro JSON encoding.
class Avro a where
    -- | Produces the schema for this type.  Schemas are represented as
    -- JSON values.  The monad is used to keep a table of already
    -- defined types, so the schema can refer to them by name.  (The
    -- concrete argument serves to specify the type, it is not actually
    -- used.)
    toSchema :: a -> MkSchema Value

    -- | Serializes a value to the binary representation.  The schema is
    -- implied, serialization to related schemas is not supported.
    toBin :: a -> Builder

    -- | Deserializzes a value from binary representation.  Right now,
    -- no attempt at schema matching is done, the schema must match the
    -- expected one exactly.
    fromBin :: Get a

    -- | Serializes a value to the JSON representation.  Note that even
    -- the JSON format needs a schema for successful deserialization,
    -- and here we support only the one implied schema.
    toAvron :: a -> Value


newtype MkSchema a = MkSchema
    { mkSchema :: (a -> H.HashMap T.Text Value -> Value) -> H.HashMap T.Text Value -> Value }

instance Functor MkSchema where fmap f m = MkSchema (\k -> mkSchema m (k . f))
instance Applicative MkSchema where pure a = MkSchema (\k -> k a)
                                    u <*> v = MkSchema (\k -> mkSchema u (\a -> mkSchema v (k . a)))
instance Monad MkSchema where return a = MkSchema (\k -> k a)
                              a >>= m = MkSchema (\k -> mkSchema a (\a' -> mkSchema (m a') k))

memoObject :: String -> [(T.Text,Value)] -> MkSchema Value
memoObject nm ps = MkSchema $ \k h ->
    let nm' = T.pack nm
        obj = object $ ("name" .= nm) : ps
    in case H.lookup nm' h of
        Nothing -> k obj $! H.insert nm' obj h
        Just obj' | obj == obj' -> k (String nm') h
                  | otherwise -> error $ "same type name, different schema: " ++ nm

runMkSchema :: MkSchema Value -> Value
runMkSchema x = mkSchema x postproc H.empty
  where
    -- Objects are fine as is.
    postproc (Object  o) _ = Object o
    -- Top level can't be a string, can it?  Need to wrap into the long form.
    postproc (String tp) _ = object [ "type" .= String tp ]
    -- Top level Array should be fine, too.
    postproc (Array a) _ = Array a
    -- reject anything else
    postproc v _ = error $ "Not allowed as toplevel schema: " ++ show v

-- instances for primitive types

-- | The Avro \"null\" type is represented as the empty tuple.
instance Avro () where
    toSchema _ = return $ String "null"
    toBin   () = mempty
    fromBin    = return ()
    toAvron () = Null

instance Avro Bool where
    toSchema _ = return $ String "boolean"
    toBin      = singleton . fromIntegral . fromEnum
    fromBin    = toEnum . fromIntegral <$> getWord8
    toAvron    = Bool

instance Avro Int where
    toSchema _ = return $ String "long"
    toBin      = encodeIntBase128
    fromBin    = decodeIntBase128
    toAvron    = Number . fromIntegral

instance Avro Int64 where
    toSchema _ = return $ String "long"
    toBin      = encodeIntBase128
    fromBin    = decodeIntBase128
    toAvron    = Number . fromIntegral

instance Avro Float where
    toSchema _ = return $ String "float"
    toBin      = putWord32le . floatToWord
    fromBin    = wordToFloat <$> getWord32le
    toAvron    = Number . fromFloatDigits

instance Avro Double where
    toSchema _ = return $ String "double"
    toBin      = putWord64le . doubleToWord
    fromBin    = wordToDouble <$> getWord64le
    toAvron    = Number . fromFloatDigits

instance Avro B.ByteString where
    toSchema _ = return $ String "bytes"
    toBin    s = encodeIntBase128 (B.length s) <> fromByteString s
    fromBin    = decodeIntBase128 >>= getByteString
    toAvron    = String . decodeLatin1

instance Avro T.Text where
    toSchema _ = return $ String "string"
    toBin      = toBin . encodeUtf8
    fromBin    = decodeUtf8 <$> fromBin
    toAvron    = String


-- Integer<->Float conversions, stolen from cereal.

{-# INLINE wordToFloat #-}
wordToFloat :: Word32 -> Float
wordToFloat x = runST (cast x)

{-# INLINE wordToDouble #-}
wordToDouble :: Word64 -> Double
wordToDouble x = runST (cast x)

{-# INLINE floatToWord #-}
floatToWord :: Float -> Word32
floatToWord x = runST (cast x)

{-# INLINE doubleToWord #-}
doubleToWord :: Double -> Word64
doubleToWord x = runST (cast x)

{-# INLINE cast #-}
cast :: ( MArray (STUArray s) b (ST s), MArray (STUArray s) a (ST s) ) => a -> ST s b
cast x = (newArray (0 :: Int, 0) x >>= castSTUArray >>= flip readArray 0)

-- | Implements Zig-Zag-Coding like in Protocol Buffers and Avro.
zig :: (Storable a, Bits a) => a -> a
zig x = (x `shiftL` 1) `xor` (x `shiftR` (8 * sizeOf x -1))

-- | Reverses Zig-Zag-Coding like in Protocol Buffers and Avro.
zag :: (Storable a, Bits a, Num a) => a -> a
zag x = negate (x .&. 1) `xor` ((x .&. complement 1) `rotateR` 1)

-- | Encodes a word of any size using a variable length "base 128"
-- encoding.
encodeWordBase128 :: (Integral a, Bits a) => a -> Builder
encodeWordBase128 x | x' == 0   = singleton (fromIntegral (x .&. 0x7f))
                    | otherwise = singleton (fromIntegral (x .&. 0x7f .|. 0x80))
                                  <> encodeWordBase128 x'
  where x' = x `shiftR` 7

decodeWordBase128 :: (Integral a, Bits a) => Get a
decodeWordBase128 = go 0 0
  where
    go acc sc = do x <- getWord8
                   let !acc' = acc .|. fromIntegral x `shiftL` sc
                   if x .&. 0x80 == 0
                        then return acc'
                        else go acc' (sc+7)

-- | Encodes an int of any size by combining the zig-zag coding with the
-- base 128 encoding.
encodeIntBase128 :: (Integral a, Bits a, Storable a) => a -> Builder
encodeIntBase128 = encodeWordBase128 . zig

-- | Decodes an int of any size by combining the zig-zag decoding with
-- the base 128 decoding.
decodeIntBase128 :: (Integral a, Bits a, Storable a) => Get a
decodeIntBase128 = zag <$> decodeWordBase128

zigInt :: Int -> Builder
zigInt = encodeIntBase128

zagInt :: Get Int
zagInt = decodeWordBase128

-- Complex Types

-- | A list becomes an Avro array
-- The chunked encoding for lists may come in handy.  How to select the
-- chunk size is not obvious, though.
instance Avro a => Avro [a] where
    toSchema as = do sa <- toSchema (head as)
                     return $ object [ "type" .= String "array", "items" .= sa ]
    toBin    [] = singleton 0
    toBin    as = toBin (length as) <> foldMap toBin as <> singleton 0
    toAvron     = Array . V.fromList . map toAvron

    -- This is not suitable for incremental processing.
    fromBin     = get_blocks []
      where
        get_blocks acc = zagInt >>= \l -> if l == 0 then return $ reverse acc
                                                    else get_block acc l >>= get_blocks
        get_block acc l = if l == 0 then return acc
                                    else fromBin >>= \a -> get_block (a:acc) (l-1)


-- | A generic vector becomes an Avro array
instance Avro a => Avro (V.Vector a) where
    toSchema as = do sa <- toSchema (V.head as)
                     return $ object [ "type" .= String "array", "items" .= sa ]
    toBin    as | V.null as = singleton 0
                | otherwise = toBin (V.length as) <> foldMap toBin as <> singleton 0
    toAvron     = Array . V.map toAvron

    -- This is not suitable for incremental processing.
    fromBin     = get_blocks []
      where
        get_blocks acc = zagInt >>= \l -> if l == 0 then return $ V.concat $ reverse acc
                                                    else get_block [] l >>=
                                                         get_blocks . (: acc) . V.fromListN l . reverse
        get_block acc l = if l == 0 then return acc
                                    else fromBin >>= \a -> get_block (a:acc) (l-1)

-- | An unboxed vector becomes an Avro array
instance (Avro a, U.Unbox a) => Avro (U.Vector a) where
    toSchema as = do sa <- toSchema (U.head as)
                     return $ object [ "type" .= String "array", "items" .= sa ]
    toBin    as | U.null as = singleton 0
                | otherwise = toBin (U.length as) <> U.foldr ((<>) . toBin) mempty as <> singleton 0
    toAvron     = Array . V.map toAvron . U.convert

    -- This is not suitable for incremental processing.
    fromBin     = get_blocks []
      where
        get_blocks acc = zagInt >>= \l -> if l == 0 then return $ U.concat $ reverse acc
                                                    else get_block [] l >>=
                                                         get_blocks . (: acc) . U.fromListN l . reverse
        get_block acc l = if l == 0 then return acc
                                    else fromBin >>= \a -> get_block (a:acc) (l-1)


-- | A map from Text becomes an Avro map.
instance Avro a => Avro (H.HashMap T.Text a) where
    toSchema   m = do sa <- toSchema (m H.! T.empty)
                      return $ object [ "type" .= String "map", "values" .= sa ]
    toBin     as | H.null as = singleton 0
                 | otherwise = toBin (H.size as) <> H.foldrWithKey (\k v b -> toBin k <> toBin v <> b) (singleton 0) as
    toAvron      = Object . H.map toAvron

    -- This is not suitable for incremental processing.
    fromBin     = get_blocks H.empty
      where
        get_blocks !acc = zagInt >>= \l -> if l == 0 then return acc
                                                     else get_block acc l >>= get_blocks
        get_block !acc l = if l == 0 then return acc
                                     else fromBin >>= \k -> fromBin >>= \v -> get_block (H.insert k v acc) (l-1)



-- * Some(!) complex types.
--
-- Enums:  Enumerated symbols.  This is generated automatically for sums
-- of empty alternatives.  Constructor names become enum symbols.

-- Records:  This is generated automatically for product types using
-- Haskell record syntax.
--
-- Unions:  For Haskell sum-of-product types using record syntax for
-- every arm, an Avro instance resolving to a union of record can be
-- generated automatically.  The constructor names become record type
-- names, their fields become record fields.

-- XXX Sometimes we build sum types containing sum types, Maybe being the
-- most obvious example.  A (Maybe a) where a itself yields a union,
-- should probably yield a union with one more alternative (the null).


deriveAvros :: [Name] -> Q [Dec]
deriveAvros = liftM concat . mapM deriveAvro

deriveAvro :: Name -> Q [Dec]
deriveAvro nm = reify nm >>= case_info
  where
    err m = fail $ "cannot derive Avro for " ++ show nm ++ ", " ++ m

    case_info (TyConI dec) = case_dec dec
    case_info            _ = err "it is not a type constructor"

    simple_cons (NormalC _ []) = True
    simple_cons _              = False

    record_cons (RecC _ _) = True
    record_cons _          = False

    case_dec (NewtypeD _cxt _name _tyvarbndrs  _con _) = err $ "don't know what to do for NewtypeD"
    case_dec (DataD    _cxt _name _tyvarbndrs cons _)
        | all simple_cons cons = mk_enum_inst [ nm1 | NormalC nm1 [] <- cons ]
        | all record_cons cons = mk_record_inst [ (nm1, vsts) | RecC nm1 vsts <- cons ]
        | otherwise            = err $ "don't know how to make an instance with these constructors"
    case_dec _ = fail $ "is not a data or newtype declaration"

    tolit = litE . StringL . nameBase
    tolitlist (x:xs) = [| T.pack $(tolit x) : $(tolitlist xs) |]
    tolitlist [    ] = [| [] |]

    -- enum instance from list of names
    mk_enum_inst :: [Name] -> Q [Dec]
    mk_enum_inst nms =
        [d| instance Avro $(conT nm) where
                toSchema _ = return $ object [ "type" .= string "enum"
                                             , "name" .= string $(tolit nm)
                                             , "symbols" .= $(tolitlist nms) ]
                toBin x = $(
                    return $ CaseE (VarE 'x)
                        [ Match (ConP nm1 [])
                                (NormalB (AppE (VarE 'zigInt)
                                               (LitE (IntegerL i)))) []
                        | (i,nm1) <- zip [0..] nms ] )

                fromBin = zagInt >>= \x -> $(
                    return $ CaseE (VarE 'x)
                        [ Match (LitP (IntegerL i))
                                (NormalB (AppE (VarE 'return)
                                               (ConE nm1))) []
                        | (i,nm1) <- zip [0..] nms ] )

                toAvron x = $(
                    return $ CaseE (VarE 'x)
                        [ Match (ConP nm1 [])
                                (NormalB (AppE (VarE 'string)
                                               (LitE (StringL (nameBase nm1))))) []
                        | nm1 <- nms ] )
        |]

    -- record instance from record-like constructors
    -- XXX maybe allow empty "normal" constructors, too
    mk_record_inst :: [ (Name, [(Name, Strict, Type)]) ] -> Q [Dec]
    mk_record_inst [(nm1,fs1)] =
        [d| instance Avro $(conT nm) where
                toSchema _ = $(mk_product_schema nm1 fs1)
                toBin      = $(to_bin_product fs1)
                fromBin    = $(from_bin_product [| return $(conE nm1) |] fs1)
                toAvron    = $(to_avron_product fs1)
        |]

    mk_record_inst arms =
        [d| instance Avro $(conT nm) where
                toSchema _ = Array . V.fromList <$> sequence
                             $( foldr (\(nm1,fs) k -> [| $(mk_product_schema nm1 fs) : $k |])
                                      [| [] |] arms )
                toBin =
                    $( do x <- newName "x"
                          LamE [VarP x] . CaseE (VarE x)
                             <$> sequence [ ($ []) . Match (RecP nm1 []) . NormalB
                                                <$> [| zigInt $(litE (IntegerL i)) <> $(to_bin_product fs) $(varE x) |]
                                          | (i,(nm1,fs)) <- zip [0..] arms ] )

                fromBin = zagInt >>=
                    $( do x <- newName "x"
                          LamE [VarP x] . CaseE (VarE x)
                            <$> sequence [ ($ []) . Match (LitP (IntegerL i)) . NormalB
                                                <$> from_bin_product [| return $(conE nm1) |] fs
                                         | (i,(nm1,fs)) <- zip [0..] arms ] )

                toAvron =
                    $( do x <- newName "x"
                          LamE [VarP x] . CaseE (VarE x)
                             <$> sequence [ ($ []) . Match (RecP nm1 []) . NormalB
                                                <$> [| object [ $(tolit nm1) .= $(to_avron_product fs) $(varE x) ] |]
                                          | (nm1,fs) <- arms ] )
        |]

    -- create schema for a product from a name and a list of fields
    mk_product_schema nm1 tps =
        [| $( fieldlist tps ) >>= \flds ->
           memoObject $( tolit nm1 )
               [ "type" .= string "record"
               , "fields" .= Array (V.fromList flds) ] |]

    fieldlist = foldr go [| return [] |]
        where
            go (nm1,_,tp) k =
                [| do sch <- toSchema $(sigE (varE 'undefined) (return tp))
                      obs <- $k
                      return $ object [ "name" .= string $(tolit nm1)
                                      , "type" .= sch ]
                             : obs |]

    -- binary encoding of records: field by field.
    to_bin_product nms =
        [| \x -> $( foldr (\(nm1,_,_) k -> [| mappend (toBin ($(varE nm1) x)) $k |] )
                          [| mempty |] nms ) |]

    from_bin_product =
        foldl (\expr (_,_,_) -> [| $expr <*> fromBin |])

    -- json encoding of records: fields in an object
    to_avron_product nms =
        [| \x -> object $(
            foldr (\(nm1,_,_) k -> [| ($(tolit nm1) .= toAvron ($(varE nm1) x)) : $k |] )
                  [| [] |] nms ) |]


data ContainerOpts = ContainerOpts { objects_per_block :: Int
                                   , filetype_label :: B.ByteString }

-- Writing a container file.  This is an 'Enumeratee', we read a list of
-- suitable types, we write a header containing the generated schema,
-- and a series of blocks with serialized data.
writeAvroContainer :: (MonadIO m, Nullable s, ListLike s a, Avro a)
                   => ContainerOpts -> Enumeratee s B.ByteString m r
writeAvroContainer ContainerOpts{..} out = do
        ma <- peekStream
        sync_marker <- liftIO $ B.pack <$> replicateM 16 randomIO

        let schema = encode . runMkSchema . toSchema . fromJust $ ma

            meta :: H.HashMap T.Text B.ByteString
            meta = H.fromList [( "avro.schema", B.concat $ BL.toChunks schema )
                              ,( "avro.codec", "null" )
                              ,( "biohazard.filetype", filetype_label )]

            hdr = fromByteString "Obj\1" <> toBin meta <> fromByteString sync_marker

        let enc_blocks = iterLoop $ \out' -> do (num,code) <- joinI $ takeStream objects_per_block $
                                                                foldStream (\(!n,c) o -> (n+1, c <> toBin o)) (0::Int,mempty)

                                                let code1 = toLazyByteString code
                                                    block = toBin num <> toBin (BL.length code1) <>
                                                            fromLazyByteString code1 <> fromByteString sync_marker
                                                lift (enumList (BL.toChunks $ toLazyByteString block) out')

        lift (enumList (BL.toChunks $ toLazyByteString hdr) out) >>= enc_blocks

-- XXX Possible codecs: null, zlib, snappy, lzma; all missing
-- XXX Should check schema on reading.

readAvroContainer :: (Monad m, ListLike s a, Avro a) => Enumeratee B.ByteString s m r
readAvroContainer out = do
        4 <- heads "Obj\1"  -- enough magic?
        meta <- iterGet (fromBin :: Get (H.HashMap T.Text B.ByteString))
        sync_marker <- iGetString 16

        flip iterLoop out $ \o -> do num <- iterGet zagInt
                                     sz <- iterGet fromBin
                                     o' <- joinI $ takeStream sz $ -- codec goes here
                                              convStream (LL.singleton `liftM` iterGet fromBin) o
                                     16 <- heads sync_marker
                                     return o'

-- | Repeatedly apply an 'Iteratee' to a value until end of stream.
-- Returns the final value.
iterLoop :: (Nullable s, Monad m) => (a -> Iteratee s m a) -> a -> Iteratee s m a
iterLoop it a = do e <- isFinished
                   if e then return a
                        else it a >>= iterLoop it


iterGet :: Monad m => Get a -> Iteratee B.ByteString m a
iterGet = go . runGetIncremental
  where
    go (Fail  _ _ err) = throwErr (iterStrExc err)
    go (Done rest _ a) = idone a (Chunk rest)
    go (Partial   dec) = liftI $ \ck -> case ck of
        Chunk s -> go (dec $ Just s)
        EOF  mx -> case dec Nothing of
            Fail  _ _ err -> throwErr (iterStrExc err)
            Partial     _ -> throwErr (iterStrExc "<partial>")
            Done rest _ a | B.null rest -> idone a (EOF mx)
                          | otherwise   -> idone a (Chunk rest)

