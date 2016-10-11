-- | Handling of Illumina BCL files.
-- We will support plain BCL, gzipped BCL and bgzf'ed BCL.  Plain BCL
-- starts with a cluster count (4 bytes, little-endian).  Base calls
-- follow with one byte per base:  bits [0..1] encode the base in the
-- order ACGT, bits 2..7 contain the quality score.
--
-- We will have to read from many files, so reading reasonably sized
-- blocks is imperative.  The typical bcl file on a MiSeq is 0.5MB, on a
-- HiSeq it's about 3MB.  We simply read them completely---this requires
-- 0.5-1GB of memory on a typical run, which shouldn't be a problem.
-- It's more if decompression is necessary, but still reasonable.

-- The BCLs come with a companion 'filter' file.  These start with three
-- header words:  zero, format version number, number of clusters.  The
-- remainder is one byte(!) per cluster, bit 0 is the filter flag.  We
-- expect one folder that contains the filter files and per-cycle
-- subfolders.

module Bio.Illumina.BCL where

import Bio.Prelude
import Bio.Util.Zlib                    ( decompressGzip )
import Data.Vector.Fusion.Util          ( Id )
import Data.Vector.Generic              ( unstream )

import qualified Data.ByteString                    as B
import qualified Data.ByteString.Lazy               as L
import qualified Data.ByteString.Lazy.Internal      as L ( ByteString(..) )
import qualified Data.ByteString.Unsafe             as B ( unsafeIndex )
import qualified Data.Vector.Fusion.Stream.Monadic  as S
import qualified Data.Vector.Fusion.Stream.Size     as S
import qualified Data.Vector.Unboxed                as U

newtype BCL  = BCL  (U.Vector Word8)
newtype Filt = Filt (U.Vector Word8)

-- | Reads a BCL file, which can be plain, or gzip'ed, or bgzf'ed.
-- We ignore the record count in the first quadword and convert
-- everything into a vector.  In case of length mismatch, we pad
-- liberally with zeroes.  The file is read and decompressed strictly.

readBCL :: FilePath -> IO BCL
readBCL fp = BCL <$> readVec fp 4

readFilt :: FilePath -> IO Filt
readFilt fp = Filt <$> readVec fp 12

readVec :: FilePath -> Int64 -> IO (U.Vector Word8)
readVec fp n = evaluate . vec_from_string . L.drop n .
                    decompressGzip . L.fromChunks . (:[]) =<< B.readFile fp

-- | Turns a lazy bytestring into a vector of words.  A straight
-- forward @fromList . toList@ would have done it, but this version
-- hopefully fuses.
vec_from_string :: L.ByteString -> U.Vector Word8
vec_from_string = unstream . S.concatMap stream_bs . stream_lbs
  where
    stream_bs :: B.ByteString -> S.Stream Id Word8
    stream_bs bs = S.Stream step 0 (S.Exact $ B.length bs)
      where
        step i | i == B.length bs = return $ S.Done
               | otherwise        = return $ S.Yield (B.unsafeIndex bs i) (i+1)

    stream_lbs :: L.ByteString -> S.Stream Id B.ByteString
    stream_lbs lbs = S.Stream step lbs S.Unknown
      where
        step  L.Empty       = return $ S.Done
        step (L.Chunk c cs) = return $ S.Yield c cs

