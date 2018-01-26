{-# LANGUAGE CPP #-}
module Bio.Util.Zlib ( decompressGzip ) where

import Prelude

import qualified Data.ByteString.Lazy            as L
import qualified Data.ByteString.Lazy.Internal   as L ( ByteString(..) )
import qualified Codec.Compression.Zlib.Internal as Z


-- | Decompresses Gzip or Bgzf and passes everything else on.  In
-- reality, it simply decompresses Gzip, and when done, looks for
-- another Gzip stream.  Since there is a small chance to attempt
-- decompression of an uncompressed stream, the original data is
-- returned in case of an error.
decompressGzip :: L.ByteString -> L.ByteString
decompressGzip s = case L.uncons s of
    Just (31, s') -> case L.uncons s' of
        Just (139,_) -> Z.foldDecompressStreamWithInput L.Chunk decompressGzip (const s)
                        (Z.decompressST Z.gzipOrZlibFormat Z.defaultDecompressParams) s
        _            -> s
    _                -> s

