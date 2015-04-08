{-# LANGUAGE TemplateHaskell, OverloadedStrings #-}
module Bio.Genocall.AvroFile where

import Bio.Base
import Bio.Bam.Pileup
import Control.Applicative
import Data.Aeson
import Data.Avro hiding ((.=))
import Data.Binary.Builder
import Data.Binary.Get
import Data.Monoid

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
    { reference_name :: T.Text
    , start_position :: Int
    , called_sites :: [ GenoCallSite ] }
  deriving (Show, Eq)

data GenoCallSite = GenoCallSite
    { snp_stats         :: CallStats
    , snp_likelihoods   :: [ Int ] -- B.ByteString?
    , indel_stats       :: CallStats
    , indel_variants    :: [ IndelVariant ]
    , indel_likelihoods :: [ Int ] } -- B.ByteString?
  deriving (Show, Eq)

deriveAvros [ ''GenoCallBlock, ''GenoCallSite, ''CallStats, ''IndelVariant ]

instance Avro V_Nuc where
    toSchema        _ = return $ object [ "type" .= String "bytes", "doc" .= String "A,C,G,T" ]
    toBin   (V_Nuc v) = encodeIntBase128 (U.length v) <> U.foldr ((<>) . singleton . unN) mempty v
    fromBin           = decodeIntBase128 >>= \l -> V_Nuc . U.fromListN l . map N . B.unpack <$> getByteString l
    toAvron (V_Nuc v) = String . T.pack . map w2c . U.toList $ U.map unN v

