{-# LANGUAGE OverloadedStrings, PatternGuards, BangPatterns, NoMonomorphismRestriction #-}
{-# OPTIONS_GHC -funbox-strict-fields #-}
module Bio.Bam.Raw (
    Block(..),
    decompressBgzfBlocks,
    decompressBgzf,
    compressBgzf,

    isBam,
    isPlainBam,
    isGzipBam,
    isBgzfBam,

    decodeBam,
    getBamRaw,
    decodeAnyBam,
    decodeAnyBamFile,
    progressPos,

    concatInputs,
    concatDefaultInputs,
    mergeInputs,
    mergeDefaultInputs,
    combineCoordinates,
    combineNames,

    encodeBam,
    encodeBamWith,
    encodeBamUncompressed,

    writeRawBamFile,
    writeRawBamHandle,
    pipeRawBamOutput,

    BamrawEnumeratee,
    BamRaw,
    bamRaw,
    virt_offset,
    raw_data
) where

import Bio.Base
import Bio.Bam.Header
import Bio.Iteratee
import Bio.Iteratee.ZLib hiding ( CompressionLevel )
import Bio.Iteratee.Bgzf
import Bio.Util ( showNum )

import Control.Applicative
import Control.Monad
import Data.Binary.Builder          ( toLazyByteString )
import Data.Binary.Put
import Data.Bits                    ( Bits, shiftL, (.|.) )
import Data.Monoid
import Data.Sequence                ( (|>) )
import System.Environment           ( getArgs )
import System.IO

import qualified Control.Monad.Catch            as C
import qualified Data.ByteString                as S
import qualified Data.ByteString.Lazy.Char8     as L
import qualified Data.ByteString.Unsafe         as S
import qualified Data.Foldable                  as F
import qualified Data.Map                       as M
import qualified Data.Sequence                  as Z

-- ^ Parsers and printers for BAM.  We employ an @Iteratee@
-- interface, and we strive to keep BAM records in their encoded form.
-- This is most compact and often fasted, since it saves the time for
-- repeated decoding and encoding, if that's not strictly needed.


type BamrawEnumeratee m b = Enumeratee' BamMeta S.ByteString [BamRaw] m b

-- | Tests if a data stream is a Bam file.
-- Recognizes plain Bam, gzipped Bam and bgzf'd Bam.  If a file is
-- recognized as Bam, a decoder (suitable Enumeratee) for it is
-- returned.  This uses 'iLookAhead' internally, so it shouldn't consume
-- anything from the stream.
isBam, isEmptyBam, isPlainBam, isBgzfBam, isGzipBam :: MonadIO m
    => Iteratee S.ByteString m (Maybe (BamrawEnumeratee m a))
isBam = firstOf [ isEmptyBam, isPlainBam, isBgzfBam, isGzipBam ]
  where
    firstOf [] = return Nothing
    firstOf (k:ks) = iLookAhead k >>= maybe (firstOf ks) (return . Just)

isEmptyBam = (\e -> if e then Just (\k -> return $ k mempty) else Nothing) `liftM` isFinished

isPlainBam = (\n -> if n == 4 then Just (joinI . decompressPlain . decodeBam) else Nothing) `liftM` heads "BAM\SOH"

-- Interesting... iLookAhead interacts badly with the parallel
-- decompression of BGZF.  (The chosen interface doesn't allow the EOF
-- signal to be passed on.)  One workaround would be to run sequential
-- BGZF decompression to check if the content is BAM, but since BGZF is
-- actually GZip in disguise, the easier workaround if to use the
-- ordinary GZip decompressor.
-- (A clean workaround would be an @Alternative@ instance for
-- @Iteratee@.)
isBgzfBam  = do b <- isBgzf
                k <- if b then joinI $ enumInflate GZip defaultDecompressParams isPlainBam else return Nothing
                return $ (\_ -> (joinI . decompressBgzfBlocks . decodeBam)) `fmap` k

isGzipBam  = do b <- isGzip
                k <- if b then joinI $ enumInflate GZip defaultDecompressParams isPlainBam else return Nothing
                return $ ((joinI . enumInflate GZip defaultDecompressParams) .) `fmap` k

-- | Checks if a file contains BAM in any of the common forms, then
-- decompresses it appropriately.  We support plain BAM, Bgzf'd BAM,
-- and Gzip'ed BAM.
--
-- The recommendation for these functions is to use @decodeAnyBam@ (or
-- @decodeAnyBamFile@) for any code that can handle @BamRaw@ input, but
-- @decodeAnyBamOrSam@ (or @decodeAnyBamOrSamFile@) for code that needs
-- @BamRec@.  That way, SAM is supported automatically, and seeking will
-- be supported if possible.
decodeAnyBam :: MonadIO m => BamrawEnumeratee m a
decodeAnyBam it = do mk <- isBam ; case mk of Just  k -> k it
                                              Nothing -> fail "this isn't BAM."

decodeAnyBamFile :: (MonadIO m, MonadMask m) => FilePath -> (BamMeta -> Iteratee [BamRaw] m a) -> m (Iteratee [BamRaw] m a)
decodeAnyBamFile fn k = enumFileRandom defaultBufSize fn (decodeAnyBam k) >>= run

concatDefaultInputs :: (MonadIO m, MonadMask m) => Enumerator' BamMeta [BamRaw] m a
concatDefaultInputs it0 = liftIO getArgs >>= \fs -> concatInputs fs it0

concatInputs :: (MonadIO m, MonadMask m) => [FilePath] -> Enumerator' BamMeta [BamRaw] m a
concatInputs [        ] = \k -> enumHandle defaultBufSize stdin (decodeAnyBam k) >>= run
concatInputs (fp0:fps0) = \k -> enum1 fp0 k >>= go fps0
  where
    enum1 "-" k1 = enumHandle defaultBufSize stdin (decodeAnyBam k1) >>= run
    enum1  fp k1 = enumFile   defaultBufSize    fp (decodeAnyBam k1) >>= run

    go [       ] = return
    go (fp1:fps) = enum1 fp1 . const >=> go fps

mergeDefaultInputs :: (MonadIO m, MonadMask m)
    => (BamMeta -> Enumeratee [BamRaw] [BamRaw] (Iteratee [BamRaw] m) a)
    -> Enumerator' BamMeta [BamRaw] m a
mergeDefaultInputs (?) it0 = liftIO getArgs >>= \fs -> mergeInputs (?) fs it0

mergeInputs :: (MonadIO m, MonadMask m)
    => (BamMeta -> Enumeratee [BamRaw] [BamRaw] (Iteratee [BamRaw] m) a)
    -> [FilePath] -> Enumerator' BamMeta [BamRaw] m a
mergeInputs  _  [        ] = \k -> enumHandle defaultBufSize stdin (decodeAnyBam k) >>= run
mergeInputs (?) (fp0:fps0) = go fp0 fps0
  where
    enum1 "-" k1 = enumHandle defaultBufSize stdin (decodeAnyBam k1) >>= run
    enum1  fp k1 = enumFile defaultBufSize fp (decodeAnyBam k1) >>= run

    go fp [       ] = enum1 fp
    go fp (fp1:fps) = mergeEnums' (go fp1 fps) (enum1 fp) (?)

{-# INLINE combineCoordinates #-}
combineCoordinates :: Monad m => BamMeta -> Enumeratee [BamRaw] [BamRaw] (Iteratee [BamRaw] m) a
combineCoordinates _ = mergeSortStreams (?)
  where u ? v = if (br_rname u, br_pos u) < (br_rname v, br_pos v) then Less else NotLess

{-# INLINE combineNames #-}
combineNames :: Monad m => BamMeta -> Enumeratee [BamRaw] [BamRaw] (Iteratee [BamRaw] m) a
combineNames _ = mergeSortStreams (?)
  where u ? v = case br_qname u `compareNames` br_qname v of LT -> Less ; _ -> NotLess



-- | Encode stuff into a BAM stream.
-- We send the encoded header and reference list to output through the
-- Bgzf compressor, then receive a list of records, which we concatenate
-- and send to output, too.
--
-- It would be nice if we were able to write an index on the side.  That
-- hasn't been designed in, yet.

encodeBam :: BamMeta -> Enumeratee [BamRaw] S.ByteString IO a
encodeBam = encodeBamWith 6 -- sensible default compression level

encodeBamUncompressed :: BamMeta -> Enumeratee [BamRaw] S.ByteString IO a
encodeBamUncompressed = encodeBamWith 0

encodeBamWith :: Int -> BamMeta -> Enumeratee [BamRaw] S.ByteString IO a
encodeBamWith lv meta = eneeBam ><> compressBgzfLv lv
  where
    eneeBam  = eneeCheckIfDone (\k -> eneeBam2 . k $ Chunk (SpecialChunk header NoChunk))
    eneeBam2 = eneeCheckIfDone (liftI . put)

    put k (EOF   mx) = idone (liftI k) $ EOF mx
    put k (Chunk rs) = eneeCheckIfDone (liftI . put) . k . Chunk $ foldr (RecordChunk . raw_data) NoChunk rs

    header = S.concat . L.toChunks $ runPut putHeader

    putHeader = do putByteString "BAM\1"
                   let hdr = toLazyByteString $ showBamMeta meta
                   putWord32le $ fromIntegral $ L.length hdr
                   putLazyByteString hdr
                   putWord32le . fromIntegral . Z.length $ meta_refs meta
                   F.mapM_ putRef $ meta_refs meta

    putRef bs = do putWord32le . fromIntegral $ S.length (sq_name bs) + 1
                   putByteString $ sq_name bs
                   putWord8 0
                   putWord32le . fromIntegral $ sq_length bs


-- | Bam record in its native encoding along with virtual address.
data BamRaw = BamRaw { virt_offset :: {-# UNPACK #-} !FileOffset
                     , raw_data :: {-# UNPACK #-} !S.ByteString }

-- | Smart constructor.  Makes sure we got a at least a full record.
bamRaw :: FileOffset -> S.ByteString -> BamRaw
bamRaw o s = if good then r else error $ "broken BAM record " ++ show (S.length s, m) ++ show m
  where
    r = BamRaw o s
    good | S.length s < 32 = False
         | otherwise       = S.length s >= sum m
    m = [ 32, br_l_read_name r, l_seq, (l_seq+1) `div` 2, n_cigar * 4 ]
    n_cigar = fromIntegral (S.unsafeIndex s 12) .|. fromIntegral (S.unsafeIndex s 13) `shiftL`  8
    l_seq = getInt s 16

-- | Accessor for raw bam.
{-# INLINE br_qname #-}
{-# DEPRECATED br_qname "use unpackBam" #-}
br_qname :: BamRaw -> Seqid
br_qname r@(BamRaw _ raw) = S.unsafeTake (br_l_read_name r) $ S.unsafeDrop 32 raw

{-# INLINE br_l_read_name #-}
{-# DEPRECATED br_l_read_name "use unpackBam" #-}
br_l_read_name :: BamRaw -> Int
br_l_read_name (BamRaw _ raw) = fromIntegral $ S.unsafeIndex raw 8 - 1

-- | Load an unaligned, little-endian int.  This is probably quite slow
-- and unnecessary on some platforms.  On i386, ix86_64 and powerPC, we
-- could cast the pointer and do a direct load.  Other may have special
-- primitives.  Worth investigating?
{-# INLINE getInt #-}
getInt :: (Num a, Bits a) => S.ByteString -> Int -> a
getInt s o = fromIntegral (S.unsafeIndex s $ o+0)             .|. fromIntegral (S.unsafeIndex s $ o+1) `shiftL`  8 .|.
             fromIntegral (S.unsafeIndex s $ o+2) `shiftL` 16 .|. fromIntegral (S.unsafeIndex s $ o+3) `shiftL` 24

{-# INLINE br_rname #-}
{-# DEPRECATED br_rname "use unpackBam" #-}
br_rname :: BamRaw -> Refseq
br_rname (BamRaw _ raw) = Refseq $ getInt raw 0

{-# INLINE br_pos #-}
{-# DEPRECATED br_pos "use unpackBam" #-}
br_pos :: BamRaw -> Int
br_pos (BamRaw _ raw) = getInt raw 4

-- | Decode a BAM stream into raw entries.  Note that the entries can be
-- unpacked using @decodeBamEntry@.  Also note that this is an
-- Enumeratee in spirit, only the @BamMeta@ and @Refs@ need to get
-- passed separately.
{-# INLINE decodeBam #-}
decodeBam :: Monad m => (BamMeta -> Iteratee [BamRaw] m a) -> Iteratee Block m (Iteratee [BamRaw] m a)
decodeBam inner = do meta <- liftBlock get_bam_header
                     refs <- liftBlock get_ref_array
                     convStream getBamRaw $ inner $! merge meta refs
  where
    get_bam_header  = do magic <- heads "BAM\SOH"
                         when (magic /= 4) $ do s <- iGetString 10
                                                fail $ "BAM signature not found: " ++ show magic ++ " " ++ show s
                         hdr_len <- endianRead4 LSB
                         joinI $ takeStream (fromIntegral hdr_len) $ parserToIteratee parseBamMeta

    get_ref_array = do nref <- endianRead4 LSB
                       foldM (\acc _ -> do
                                   nm <- endianRead4 LSB >>= iGetString . fromIntegral
                                   ln <- endianRead4 LSB
                                   return $! acc |> BamSQ (S.init nm) (fromIntegral ln) []
                             ) Z.empty $ [1..nref]

    -- Need to merge information from header into actual reference list.
    -- The latter is the authoritative source for the *order* of the
    -- sequences, so leftovers from the header are discarded.  Merging
    -- is by name.  So we merge information from the header into the
    -- list, then replace the header information.
    merge meta refs =
        let tbl = M.fromList [ (sq_name sq, sq) | sq <- F.toList (meta_refs meta) ]
        in meta { meta_refs = fmap (\s -> maybe s (merge' s) (M.lookup (sq_name s) tbl)) refs }

    merge' l r | sq_length l == sq_length r = l { sq_other_shit = sq_other_shit l ++ sq_other_shit r }
               | otherwise                  = l -- contradiction in header, but we'll just ignore it


{-# INLINE getBamRaw #-}
getBamRaw :: Monad m => Iteratee Block m [BamRaw]
getBamRaw = do off <- getOffset
               raw <- liftBlock $ do
                        bsize <- endianRead4 LSB
                        when (bsize < 32) $ fail "short BAM record"
                        iGetString (fromIntegral bsize)
               return [bamRaw off raw]

writeRawBamFile :: FilePath -> BamMeta -> Iteratee [BamRaw] IO ()
writeRawBamFile fp meta =
    C.bracket (liftIO $ openBinaryFile fp WriteMode)
              (liftIO . hClose)
              (flip writeRawBamHandle meta)

pipeRawBamOutput :: BamMeta -> Iteratee [BamRaw] IO ()
pipeRawBamOutput meta = encodeBamUncompressed meta =$ mapChunksM_ (liftIO . S.hPut stdout)

writeRawBamHandle :: Handle -> BamMeta -> Iteratee [BamRaw] IO ()
writeRawBamHandle hdl meta = encodeBam meta =$ mapChunksM_ (liftIO . S.hPut hdl)

-- | A simple progress indicator that prints sequence id and position.
progressPos :: MonadIO m => String -> (String -> IO ()) -> Refs -> Enumeratee [BamRaw] [BamRaw] m a
progressPos msg put refs = eneeCheckIfDonePass (icont . go 0)
  where
    go !_ k (EOF         mx) = idone (liftI k) (EOF mx)
    go !n k (Chunk    [   ]) = liftI $ go n k
    go !n k (Chunk as@(a:_)) = do let !n' = n + length as
                                      nm = unpackSeqid (sq_name (getRef refs (br_rname a))) ++ ":"
                                  when (n `div` 65536 /= n' `div` 65536) $ liftIO $ put $
                                        "\27[K" ++ msg ++ nm ++ showNum (br_pos a) ++ "\r"
                                  eneeCheckIfDonePass (icont . go n') . k $ Chunk as

