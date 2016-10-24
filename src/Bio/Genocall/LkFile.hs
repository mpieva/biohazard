{-# LANGUAGE DeriveGeneric, StandaloneDeriving, CPP #-}
#if __GLASGOW_HASKELL__ > 710
{-# OPTIONS_GHC -Wno-orphans #-}
#else
{-# OPTIONS_GHC -fno-warn-orphans #-}
#endif

module Bio.Genocall.LkFile where

import Bio.Bam.Header
import Bio.Bam.Pileup
import Bio.Prelude
import Data.Binary.Serialise.CBOR
import Data.MiniFloat

import qualified Data.Vector.Unboxed            as U

-- ^ File format for genotype calls.
--
-- We'll write a "CBOR sequence file".  We start with a header, end with
-- a footer.  Records will not needlessly repeat coordinates, so we
-- write blocks consisting of a block header and the stuff inside.
-- Sooo... to turn the file into a list of a homogenous type, we need a
-- variant.

data GenoFileRec
        = GenoFileHeader GenoHeader
        | GenoFileBlock  GenoCallBlock
        | GenoFileSite   GenoCallSite
        | GenoFileFooter GenoFooter
  deriving (Show, Eq, Generic)


data GenoHeader = GenoHeader
        { header_version :: Int
        , header_refs :: Refs }
  deriving (Show, Eq, Generic)

data GenoFooter = GenoFooter
  deriving (Show, Eq, Generic)

-- | To output a container file, we need to convert calls into a stream of
-- sensible objects.  To cut down on redundancy, the object will have a
-- header that names the reference sequence and the start, followed by
-- calls.  The calls themselves have contiguous coordinates, we start a
-- new block if we have to skip; we also start a new block when we feel
-- the current one is getting too large.

data GenoCallBlock = GenoCallBlock
        { reference_name :: {-# UNPACK #-} !Refseq
        , start_position :: {-# UNPACK #-} !Int }
  deriving (Show, Eq, Generic)

data GenoCallSite = GenoCallSite
        { snp_stats         :: {-# UNPACK #-} !CallStats
        -- snp likelihoods appear in the same order as in VCF, the reference
        -- allele goes first if it is A, C, G or T.  Else A goes first---not
        -- my problem how to express that in VCF.
        , snp_likelihoods   :: {-# UNPACK #-} !(U.Vector Mini) -- Bytes?
        , ref_allele        :: {-# UNPACK #-} !Nucleotides
        , indel_stats       :: {-# UNPACK #-} !CallStats
        , indel_variants    ::                [ IndelVariant ]
        , indel_likelihoods :: {-# UNPACK #-} !(U.Vector Mini) } -- Bytes?
  deriving (Show, Eq, Generic)

instance Serialise GenoFileRec
instance Serialise GenoHeader
instance Serialise GenoFooter
instance Serialise GenoCallBlock
instance Serialise GenoCallSite
instance Serialise IndelVariant
instance Serialise CallStats


instance Serialise Refseq where
    encode = encode . unRefseq
    decode = Refseq <$> decode

instance Serialise Nucleotides where
    encode = encode . unNs
    decode = Ns <$> decode

instance Serialise Nucleotide where
    encode = encode . unN
    decode = N <$> decode

instance Serialise BamSQ where
    encode sq = encode (sq_name sq, sq_length sq)
    decode = (\(n,l) -> BamSQ n l []) <$> decode

instance Serialise Mini where
    encode = encode . unMini
    decode = Mini <$> decode

instance Serialise V_Nucs where
    encode (V_Nucs v) = encode v
    decode = V_Nucs <$> decode

instance Serialise V_Nuc where
    encode (V_Nuc v) = encode v
    decode = V_Nuc <$> decode

-- | Storing likelihoods:  we take the natural logarithm (GL values are
-- already in a log scale) and convert to minifloat 0.4.4
-- representation.  Range and precision should be plenty.
compact_likelihoods :: U.Vector Prob -> U.Vector Mini -- Bytes?
compact_likelihoods = U.map $ float2mini . negate . unPr
-- compact_likelihoods = map fromIntegral {- B.pack -} . U.toList . U.map (float2mini . negate . unPr)

