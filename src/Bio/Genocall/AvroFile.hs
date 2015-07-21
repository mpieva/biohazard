{-# LANGUAGE TemplateHaskell, OverloadedStrings #-}
module Bio.Genocall.AvroFile where

import Bio.Base
import Bio.Bam.Pileup
import Control.Applicative
import Data.Aeson
import Data.Avro hiding ( (.=) )
import Data.Binary.Builder
import Data.Binary.Get
import Data.List ( intersperse )
import Data.Monoid
import Data.MiniFloat

import qualified Data.ByteString                as B
import qualified Data.Text                      as T
import qualified Data.Vector.Unboxed            as U

-- ^ File format for genotype calls.

-- | To output a container file, we need to convert calls into a stream of
-- sensible objects.  To cut down on redundancy, the object will have a
-- header that names the reference sequence and the start, followed by
-- calls.  The calls themselves have contiguous coordinates, we start a
-- new block if we have to skip; we also start a new block when we feel
-- the current one is getting too large.

data GenoCallBlock = GenoCallBlock
    { reference_name :: {-# UNPACK #-} !T.Text
    , start_position :: {-# UNPACK #-} !Int
    , called_sites :: [ GenoCallSite ] }
  deriving (Show, Eq)

data GenoCallSite = GenoCallSite
    { snp_stats         :: {-# UNPACK #-} !CallStats
    , snp_likelihoods   :: {-# UNPACK #-} !(U.Vector Mini) -- B.ByteString?
    , ref_allele        :: {-# UNPACK #-} !Nucleotides
    , indel_stats       :: {-# UNPACK #-} !CallStats
    , indel_variants    :: [ IndelVariant ]
    , indel_likelihoods :: {-# UNPACK #-} !(U.Vector Mini) } -- B.ByteString?
  deriving (Show, Eq)

-- | Storing likelihoods:  we take the natural logarithm (GL values are
-- already in a log scale) and convert to minifloat 0.4.4
-- representation.  Range and precision should be plenty.
compact_likelihoods :: U.Vector (Prob Double) -> U.Vector Mini -- B.ByteString
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

