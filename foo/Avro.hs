{-# LANGUAGE OverloadedStrings, FlexibleInstances, TemplateHaskell #-}
module Avro where

import Control.Applicative
import Control.Monad
import Data.Aeson
import Data.Bits
import Data.ByteString.Builder
import Data.Foldable ( foldMap )
import Data.Monoid
import Data.Scientific
import Data.Text.Encoding
import Foreign.Storable ( Storable, sizeOf )
import Language.Haskell.TH

import qualified Data.ByteString as B
import qualified Data.ByteString.Lazy as BL
import qualified Data.HashMap.Strict as H
import qualified Data.Text as T
import qualified Data.Text.Lazy as TL
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U


-- | This is the class of types we can embed into the Avro
-- infrastructure.  Right now, we can derive a schema, toBin     to
-- the Avro binary format, and toBin     to the Avro JSON encoding.
-- XXX Deserialization is out of scope right now.
class Avro a where
    -- | Produces the schema for this type.  Schemas are represented as
    -- JSON values.  The monad is used to keep a table of already
    -- defined types, so the schema can refer to them by name.  (The
    -- concrete argument serves to specify the type, it is not actually
    -- used.)
    toSchema :: a -> MkSchema Value

    -- | Serializes a value to the binary representation.  The schema is
    -- implied, serialization to related schemas is not supported.
    toBin     :: a -> Builder

    -- | Serializes a value to the JSON representation.  Note that even
    -- the JSON format needs a schema for successful deserialization,
    -- and here we support only the one implied schema.
    toAvron  :: a -> Value


newtype MkSchema a = MkSchema a -- XXX

instance Functor MkSchema where fmap f (MkSchema a) = MkSchema (f a)
instance Monad MkSchema where return = MkSchema
                              MkSchema a >>= k = k a

memoObject :: T.Text -> [(T.Text,Value)] -> MkSchema Value
memoObject nm ps = return . object $ ("name" .= nm) : ps

-- instances for primitive types

-- | The Avro \"null\" type is represented as the empty tuple.
instance Avro () where
    toSchema _ = return $ String "null"
    toBin   () = mempty
    toAvron () = Null

instance Avro Bool where
    toSchema _ = return $ String "boolean"
    toBin      = word8 . fromIntegral . fromEnum
    toAvron    = Bool

instance Avro Int where
    toSchema _ = return $ String "long"
    toBin      = encodeIntBase128
    toAvron    = Number . fromIntegral

instance Avro Float where
    toSchema _ = return $ String "float"
    toBin      = floatLE
    toAvron    = Number . fromFloatDigits

instance Avro Double where
    toSchema _ = return $ String "double"
    toBin      = doubleLE
    toAvron    = Number . fromFloatDigits

instance Avro B.ByteString where
    toSchema _ = return $ String "bytes"
    toBin    s = encodeIntBase128 (B.length s) <> byteString s
    toAvron    = String . decodeLatin1

instance Avro T.Text where
    toSchema _ = return $ String "string"
    toBin      = toBin . encodeUtf8
    toAvron    = String




-- | Implements Zig-Zag-Coding like in Protocol Buffers and Avro.
zig :: (Storable a, Bits a) => a -> a
zig x = (x `shiftL` 1) `xor` (x `shiftR` (8 * sizeOf x -1))

-- | Reverses Zig-Zag-Coding like in Protocol Buffers and Avro.
zag :: (Storable a, Bits a, Num a) => a -> a
zag x = negate (x .&. 1) `xor` ((x .&. complement 1) `rotateR` 1)

-- | Encodes a word of any size using a variable length "base 128"
-- encoding.
encodeWordBase128 :: (Integral a, Bits a) => a -> Builder
encodeWordBase128 x | x' == 0   = word8 (fromIntegral (x .&. 0x7f))
                    | otherwise = word8 (fromIntegral (x .&. 0x7f .|. 0x80))
                                  <> encodeWordBase128 x'
  where x' = x `shiftR` 7

-- | Encodes an int of any size by combining the zig-zag coding with the
-- base 128 encoding.
encodeIntBase128 :: (Integral a, Bits a, Storable a) => a -> Builder
encodeIntBase128 = encodeWordBase128 . zig

zigInt :: Int -> Builder
zigInt = encodeIntBase128

-- Complex Types

-- | A list becomes an Avro array
-- The chunked encoding for lists may come in handy.  How to select the
-- chunk size is not obvious, though.
instance Avro a => Avro [a] where
    toSchema as = do sa <- toSchema (head as)
                     return $ object [ "type" .= String "array", "items" .= sa ]
    toBin    as = toBin (length as) <> foldMap toBin as <> word8 0
    toAvron     = Array . V.fromList . map toAvron

-- | A generic vector becomes an Avro array
instance Avro a => Avro (V.Vector a) where
    toSchema as = do sa <- toSchema (V.head as)
                     return $ object [ "type" .= String "array", "items" .= sa ]
    toBin    as = toBin (V.length as) <> foldMap toBin as <> word8 0
    toAvron     = Array . V.map toAvron

-- | An unboxed vector becomes an Avro array
instance (Avro a, U.Unbox a) => Avro (U.Vector a) where
    toSchema as = do sa <- toSchema (U.head as)
                     return $ object [ "type" .= String "array", "items" .= sa ]
    toBin    as = toBin (U.length as) <> U.foldr ((<>) . toBin) mempty as <> word8 0
    toAvron     = Array . V.map toAvron . U.convert

-- | A map from Text becomes an Avro map.
instance Avro a => Avro (H.HashMap T.Text a) where
    toSchema   m = do sa <- toSchema (m H.! T.empty)
                      return $ object [ "type" .= String "map", "values" .= sa ]
    toBin     as = toBin (H.size as) <> H.foldrWithKey (\k v b -> toBin k <> toBin v <> b) (word8 0) as
    toAvron      = Object . H.map toAvron


-- Enums
-- Enumerated symbols.  Sums of empty alternatives would fit.
-- Can we do this generically?  I'm too retarded to do it.

-- Unions
-- Unions are unions of differently named types.  Note that unions
-- themselves are anonymous.

-- Fixed
-- Fixed size, uninterpreted things.  Could come in handy for, say,
-- exactly 10 genotype likelihoods.

-- Records
-- Can we do this generically?  I'm too retarded to do it.


-- We could encode a typical sum-of products as a Union of Records.  The
-- constructor names become record type names, their fields become
-- record fields.  Could be done uding TH?

-- The degenerate case of a single data constructor can be treated
-- specially, it doesn't need to be a union.

-- Sometimes we build sum types containing sum types, Maybe being the
-- most obvious example.  A (Maybe a) where a itself yields a union,
-- should probably yield a union with one more alternative (the null).

-- deriving stuff for stuff...
--
-- The Dec inside must be a DataD or a NewtypeD.  Ignore Cxt, ignore
-- Name, ignore(?) TyVarBndrs, ignore deriving Names.  Recurse into
-- Cons.
--
-- If all Cons are RecC: generate an Avro record type from the Name,
-- with one field per VarStrictType.  Name becomes field name, Type is
-- ignored (we recurse).  We'll need to add an (Avro t) context if t
-- contains a variable.
--
-- If all Cons are NormalC with no arguments, create an Enum from their
-- names.  


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
        | all simple_cons cons = mk_enum_inst [ nm | NormalC nm [] <- cons ]
        | all record_cons cons = mk_record_inst [ (nm, vsts) | RecC nm vsts <- cons ]
        | otherwise            = err $ "don't know how to make an instance with these constructors"
    case_dec _ = fail $ "is not a data or newtype declaration"

    tolit = litE . StringL . nameBase
    tolitlist (x:xs) = [| T.pack $(tolit x) : $(tolitlist xs) |]
    tolitlist [    ] = [| [] |]

    -- enum instance from list of names
    mk_enum_inst :: [Name] -> Q [Dec]
    mk_enum_inst nms =
        [d| instance Avro $(conT nm) where
                toSchema _ = return $ object [ T.pack "type" .= T.pack "enum"
                                             , T.pack "name" .= T.pack $(tolit nm)
                                             , T.pack "symbols" .= $(tolitlist nms) ]
                toBin x = $( 
                    return $ CaseE (VarE 'x) 
                        [ Match (ConP nm1 [])
                                (NormalB (AppE (VarE 'zigInt)
                                               (LitE (IntegerL i)))) []
                        | (i,nm1) <- zip [0..] nms ] )
                toAvron x = $(
                    return $ CaseE (VarE 'x) 
                        [ Match (ConP nm1 [])
                                (NormalB (AppE (ConE 'String)
                                               (AppE (VarE 'T.pack)
                                                     (LitE (StringL (nameBase nm1)))))) []
                        | (i,nm1) <- zip [0..] nms ] )
        |]

    -- record instance from record-like constructors
    -- XXX maybe allow empty "normal" constructors, too
    -- XXX if there is only one alternative, the top-level union should not be generated
    -- each Name becomes one alternative, each VarStrictType becomes a field.
    mk_record_inst :: [ (Name, [(Name, Strict, Type)]) ] -> Q [Dec]
    mk_record_inst [(nm1,fs1)] =
        [d| instance Avro $(conT nm) where
                toSchema _ = $(mk_product_schema nm1 fs1)
                toBin      = $(to_bin_product fs1)
                toAvron    = $(to_avron_product fs1)
        |]

    mk_record_inst arms =
        [d| instance Avro $(conT nm) where
                {- toSchema _ = return . Array $ V.fromList 
                             $( mk_record_schema object [ T.pack "type" .= "union"
                                             , T.pack "name" .= $(tolit nm) ]
                                             -- , T.pack "symbols" .= $(tolitlist nms) ] -}
                toBin = 
                    $( do x <- newName "x"
                          LamE [VarP x] . CaseE (VarE x) 
                             <$> sequence [ ($ []) . Match (RecP nm []) . NormalB
                                                <$> [| zigInt $(litE (IntegerL i)) <> $(to_bin_product fs) $(varE x) |]
                                          | (i,(nm,fs)) <- zip [0..] arms ] )
                toAvron = 
                    $( do x <- newName "x"
                          LamE [VarP x] . CaseE (VarE x) 
                             <$> sequence [ ($ []) . Match (RecP nm []) . NormalB
                                                <$> [| object [ T.pack $(tolit nm) .= $(to_avron_product fs) $(varE x) ] |]
                                          | (nm,fs) <- arms ] )
        |]

    -- create schema for a product from a name and a list of fields
    mk_product_schema nm tps =
        [| $( fieldlist tps ) >>= \flds ->
           memoObject $( tolit nm )
               [ "type" .= String "record"
               , "fields" .= Array (V.fromList flds) ] |]

    fieldlist = foldr go [| return [] |]
        where
            go (nm,_,tp) k = 
                [| do sch <- toSchema $(sigE (varE 'undefined) (return tp))
                      obs <- $k
                      return $ object [ "name" .= T.pack $(tolit nm)
                                      , "type" .= sch ]
                             : obs |]

    -- binary encoding of records: field by field.
    to_bin_product nms = 
        [| \x -> $( foldr (\(nm,_,_) k -> [| mappend (toBin ($(varE nm) x)) $k |] )
                          [| mempty |] nms ) |]

    -- json encoding of records: fields in an object
    to_avron_product nms =
        [| \x -> object $( 
            foldr (\(nm,_,_) k -> [| (T.pack $(tolit nm) .= toAvron ($(varE nm) x)) : $k |] )
                  [| [] |] nms ) |]
