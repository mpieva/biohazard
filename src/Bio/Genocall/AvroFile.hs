{-# LANGUAGE TemplateHaskell, OverloadedStrings, PatternGuards #-}
module Bio.Genocall.AvroFile where

import Bio.Bam.Header
import Bio.Bam.Pileup
import Bio.Prelude
import Data.Aeson
import Data.Avro
import Data.Binary.Builder
import Data.Binary.Get
import Data.MiniFloat
import Data.Scientific ( toBoundedInteger )
import Data.Text.Encoding ( encodeUtf8 )

import qualified Data.ByteString                as B
import qualified Data.HashMap.Strict            as H
import qualified Data.Text                      as T
import qualified Data.Vector                    as V
import qualified Data.Vector.Unboxed            as U
import qualified Data.Sequence                  as Z

-- ^ File format for genotype calls.

-- | To output a container file, we need to convert calls into a stream of
-- sensible objects.  To cut down on redundancy, the object will have a
-- header that names the reference sequence and the start, followed by
-- calls.  The calls themselves have contiguous coordinates, we start a
-- new block if we have to skip; we also start a new block when we feel
-- the current one is getting too large.

data GenoCallBlock = GenoCallBlock
    { reference_name :: {-# UNPACK #-} !Refseq
    , start_position :: {-# UNPACK #-} !Int
    , called_sites :: [ GenoCallSite ] }
  deriving (Show, Eq)

data GenoCallSite = GenoCallSite
    { snp_stats         :: {-# UNPACK #-} !CallStats
    -- snp likelihoods appear in the same order as in VCF, the reference
    -- allele goes first if it is A, C, G or T.  Else A goes first---not
    -- my problem how to express that in VCF.
    , snp_likelihoods   :: {-# UNPACK #-} !(U.Vector Mini) -- B.ByteString?
    , ref_allele        :: {-# UNPACK #-} !Nucleotides
    , indel_stats       :: {-# UNPACK #-} !CallStats
    , indel_variants    :: [ IndelVariant ]
    , indel_likelihoods :: {-# UNPACK #-} !(U.Vector Mini) } -- B.ByteString?
  deriving (Show, Eq)

-- | Storing likelihoods:  we take the natural logarithm (GL values are
-- already in a log scale) and convert to minifloat 0.4.4
-- representation.  Range and precision should be plenty.
compact_likelihoods :: U.Vector Prob -> U.Vector Mini -- B.ByteString
compact_likelihoods = U.map $ float2mini . negate . unPr
-- compact_likelihoods = map fromIntegral {- B.pack -} . U.toList . U.map (float2mini . negate . unPr)


deriveAvros [ ''GenoCallBlock, ''GenoCallSite, ''CallStats, ''IndelVariant ]

instance Avro V_Nuc where
    toSchema        _ = return $ object [ "type" .= String "bytes", "doc" .= String doc ]
      where doc = T.pack $ intersperse ',' $ show $ [minBound .. maxBound :: Nucleotide]
    toBin   (V_Nuc v) = encodeIntBase128 (U.length v) <> U.foldr ((<>) . singleton . unN) mempty v
    fromBin           = decodeIntBase128 >>= \l -> V_Nuc . U.fromListN l . map N . B.unpack <$> getByteString l
    toAvron (V_Nuc v) = String . T.pack . map w2c . U.toList $ U.map unN v

instance Avro V_Nucs where
    toSchema         _ = return $ object [ "type" .= String "bytes", "doc" .= String doc ]
      where doc = T.pack $ intersperse ',' $ show $ [minBound .. maxBound :: Nucleotides]
    toBin   (V_Nucs v) = encodeIntBase128 (U.length v) <> U.foldr ((<>) . singleton . unNs) mempty v
    fromBin            = decodeIntBase128 >>= \l -> V_Nucs . U.fromListN l . map Ns . B.unpack <$> getByteString l
    toAvron (V_Nucs v) = String . T.pack . map w2c . U.toList $ U.map unNs v

instance Avro Nucleotides where
    toSchema _ = return $ String "int"
    toBin      = encodeIntBase128 . unNs
    fromBin    = Ns <$> decodeIntBase128
    toAvron    = Number . fromIntegral . unNs

instance Avro Mini where
    toSchema _ = return $ String "int"
    toBin      = encodeIntBase128 . unMini
    fromBin    = Mini <$> decodeIntBase128
    toAvron    = Number . fromIntegral . unMini

-- | We encode the Refseq as an Avro enum, which serves as a kind of
-- symbol table.  To make this work, the environment of the 'MkSchema'
-- monad has to be prepopulated with a suitable schema.
instance Avro Refseq where
    toSchema _ = getNamedSchema "Refseq"
    toBin      = encodeIntBase128 . unRefseq
    fromBin    = Refseq <$> decodeIntBase128

    -- This is cheating, we should use the enum names, but they are not
    -- available.  Doesn't matter, this is mostly for debugging anyway.
    toAvron    = Number . fromIntegral . unRefseq


-- | Reconstructs the list of reference sequences from Avro metadata.
-- If a type named @Refseq@ is defined in the schema and is an enum, it
-- defines the symbol table, otherwise an empty list is returned.  If
-- @biohazard.refseq_length@ exists, and is an array, it's elements are
-- interpreted as the lengths in order, otherwise the lengths are set to
-- zero.
getRefseqs :: AvroMeta -> Refs
getRefseqs meta
    | Object o <- findSchema "Refseq" meta
    , Just (String "enum") <- H.lookup "type" o
    , Just (Array    syms) <- H.lookup "symbols" o
            = Z.fromList [ BamSQ (encodeUtf8 nm) ln [] | (String nm, ln) <- V.toList syms `zip` lengths ]
    | otherwise = Z.empty
  where
    lengths = case decodeStrict =<< H.lookup "biohazard.refseq_length" meta of
        Just (Array lns) -> [ case l of Number n -> maybe 0 id $ toBoundedInteger n ; _ -> 0 | l <- V.toList lns ]
        _                -> repeat 0

