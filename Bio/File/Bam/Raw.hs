{-# LANGUAGE OverloadedStrings, PatternGuards, BangPatterns, NoMonomorphismRestriction #-}
module Bio.File.Bam.Raw (
    module Bio.Base,

    Block(..),
    decompressBgzf,
    compressBgzf,

    isBam,
    isPlainBam,
    isGzipBam,
    isBgzfBam,

    decodeBam,
    decodeBamSequence,
    decodeAnyBam,
    decodeAnyBamFile,

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
    raw_data,
    br_qname,
    br_l_read_name,
    br_l_seq,
    br_n_cigar_op,
    br_aln_length,
    br_flag,

    br_isPaired,
    br_isProperlyPaired,
    br_isUnmapped,
    br_isMateUnmapped,
    br_isReversed,
    br_isMateReversed,
    br_isFirstMate,
    br_isSecondMate,
    br_isAuxillary,
    br_isFailsQC,
    br_isDuplicate,

    br_rname,
    br_pos,
    br_mrnm,
    br_mpos,
    br_isize,

    Refs,
    noRefs,
    getRef,

    Refseq(..),
    invalidRefseq,
    isValidRefseq,
    invalidPos,
    isValidPos,
    unknownMapq,

    BamIndex,
    readBamIndex,
    readBamIndex',

    Mutator,
    mutateBamRaw,
    setFlag,
    setMrnm,
    setMpos,
    setIsize
) where

import Bio.Base
import Bio.File.Bam.Header
import Bio.File.Bgzf
import Bio.Iteratee
import Bio.Iteratee.ZLib hiding ( CompressionLevel )

import Control.Monad
import Data.Array.IO                ( IOUArray, newArray_, writeArray )
import Data.Array.Unboxed
import Data.Array.Unsafe            ( unsafeFreeze )
import Data.Attoparsec.Iteratee
import Data.Binary.Put
import Data.Bits                    ( Bits, shiftL, shiftR, (.&.), (.|.), testBit )
import Data.Int                     ( Int64 )
import Data.Monoid
import Data.Sequence                ( (|>) )
import Foreign.C.String             ( CString )
import Foreign.Storable             ( pokeElemOff )
import System.Environment           ( getArgs )
import System.IO
import System.IO.Unsafe

import qualified Control.Monad.CatchIO          as CIO
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


-- | The invalid position.
-- Bam uses this value to encode a missing position.
invalidPos :: Int
invalidPos = 0xFFFFFFFF

-- | Tests whether a position is valid.
-- Returns true unless the the argument equals @invalidPos@.
isValidPos :: Int -> Bool
isValidPos = (/=) invalidPos

unknownMapq :: Int
unknownMapq = 255

type BamrawEnumeratee m b = Enumeratee' BamMeta S.ByteString [BamRaw] m b

-- | Tests if a data stream is a Bam file.
-- Recognizes plain Bam, gzipped Bam and bgzf'd Bam.  If a file is
-- recognized as Bam, a decoder (suitable Enumeratee) for it is returned.
isBam, isEmptyBam, isPlainBam, isBgzfBam, isGzipBam :: MonadIO m
    => Iteratee S.ByteString m (Maybe (BamrawEnumeratee m a))
isBam = firstOf [ isEmptyBam, isPlainBam, isBgzfBam, isGzipBam ]
  where
    firstOf [] = return Nothing
    firstOf (k:ks) = k >>= maybe (firstOf ks) (return . Just)

isEmptyBam = (\e -> if e then Just (\k -> return $ k mempty) else Nothing) `liftM` isFinished

isPlainBam = (\n -> if n == 4 then Just (joinI . decompressPlain . decodeBam) else Nothing) 
             `liftM` i'lookAhead (heads "BAM\SOH")

-- Interesting... i'lookAhead interacts badly with the parallel
-- decompression of BGZF.  (The chosen interface doesn't allow the EOF
-- signal to be passed on.)  One workaround would be to run sequential
-- BGZF decompression to check if the content is BAM, but since BGZF is
-- actually GZip in disguise, the easier workaround if to use the
-- ordinary GZip decompressor.
-- (A clean workaround would be an @Alternative@ instance for
-- @Iteratee@.)
isBgzfBam  = do b <- isBgzf
                k <- if b then i'lookAhead $ joinI $ enumInflate GZip defaultDecompressParams isPlainBam
                          else return Nothing
                return $ (\_ -> (joinI . decompressBgzf . decodeBam)) `fmap` k

isGzipBam  = do b <- isGzip
                k <- if b then i'lookAhead $ joinI $ enumInflate GZip defaultDecompressParams isPlainBam
                          else return Nothing
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

decodeAnyBamFile :: MonadCatchIO m => FilePath -> (BamMeta -> Iteratee [BamRaw] m a) -> m (Iteratee [BamRaw] m a)
decodeAnyBamFile fn k = enumFileRandom defaultBufSize fn (decodeAnyBam k) >>= run

-- | Seek to a given sequence in a Bam file, read those records.  This
-- requires an appropriate index (read separately), and the file must
-- have been opened in such a way as to allow seeking.  Enumerates over
-- the @BamRaw@ records of the correct sequence only, doesn't enumerate
-- at all if the sequence isn't found.

decodeBamSequence :: Monad m => BamIndex -> Refseq -> Enumeratee Block [BamRaw] m a
decodeBamSequence idx refseq iter = case idx ! refseq of
        _ | not (bounds idx `inRange` refseq) -> return iter
        0                                     -> return iter
        virtoff -> do seek $ fromIntegral virtoff
                      (decodeBamLoop ><> breakE wrong_ref) iter
  where
    wrong_ref br = let a = fromIntegral $ raw_data br `S.index` 0
                       b = fromIntegral $ raw_data br `S.index` 1
                       c = fromIntegral $ raw_data br `S.index` 2
                       d = fromIntegral $ raw_data br `S.index` 3
                       r = a `shiftL`  0 .|.  b `shiftL`  8 .|. 
                           c `shiftL` 16 .|.  d `shiftL` 24
                   in r /= unRefseq refseq

mergeDefaultInputs :: MonadCatchIO m
    => (BamMeta -> Enumeratee [BamRaw] [BamRaw] (Iteratee [BamRaw] m) a)
    -> Enumerator' BamMeta [BamRaw] m a
mergeDefaultInputs (?) it0 = liftIO getArgs >>= \fs -> mergeInputs (?) fs it0

mergeInputs :: MonadCatchIO m
    => (BamMeta -> Enumeratee [BamRaw] [BamRaw] (Iteratee [BamRaw] m) a)
    -> [FilePath] -> Enumerator' BamMeta [BamRaw] m a
mergeInputs  _  [        ] = \k -> enumHandle defaultBufSize stdin (decodeAnyBam k) >>= run
mergeInputs (?) (fp0:fps0) = go fp0 fps0
  where
    enum1 "-" k1 = enumHandle defaultBufSize stdin (decodeAnyBam k1) >>= run
    enum1  fp k1 = enumFile defaultBufSize fp (decodeAnyBam k1) >>= run

    go fp [       ] = enum1 fp
    go fp (fp1:fps) = mergeEnums' (go fp1 fps) (enum1 fp) (?)

combineCoordinates :: Monad m => BamMeta -> Enumeratee [BamRaw] [BamRaw] (Iteratee [BamRaw] m) a
combineCoordinates _ = mergeSortStreams (?)
  where u ? v = if (br_rname u, br_pos u) < (br_rname v, br_pos v) then Less else NotLess

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

encodeBam :: MonadIO m => BamMeta -> Enumeratee [BamRaw] S.ByteString m a
encodeBam = encodeBamWith 6 -- sensible default compression level

encodeBamUncompressed :: MonadIO m => BamMeta -> Enumeratee [BamRaw] S.ByteString m a
encodeBamUncompressed = encodeBamWith 0

encodeBamWith :: MonadIO m => Int -> BamMeta -> Enumeratee [BamRaw] S.ByteString m a
encodeBamWith lv meta = eneeBam ><> compressBgzf lv
  where
    eneeBam = eneeCheckIfDone (\k -> eneeBam2 . k $ Chunk header)
    eneeBam2 = eneeCheckIfDone (\k -> eneeBam3 . k $ Chunk S.empty)
    eneeBam3 = eneeCheckIfDone (liftI . put)

    put k (EOF                mx) = idone (liftI k) $ EOF mx
    put k (Chunk [             ]) = liftI $ put k
    put k (Chunk (BamRaw _ r:rs)) = eneeCheckIfDone (\k' -> put k' (Chunk rs)) . k $ Chunk r'
      where
        l  = S.length r
        r' = S.cons (fromIntegral (l `shiftR`  0 .&. 0xff)) $
             S.cons (fromIntegral (l `shiftR`  8 .&. 0xff)) $
             S.cons (fromIntegral (l `shiftR` 16 .&. 0xff)) $
             S.cons (fromIntegral (l `shiftR` 24 .&. 0xff)) r
    
    header = S.concat . L.toChunks $ runPut putHeader

    putHeader = do putByteString "BAM\1"
                   let hdr = showBamMeta meta L.empty
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
bamRaw o s = if good then r else error $ "broken BAM record " ++ show (S.length s, m) ++
    show [ 32, br_l_read_name r, br_l_seq r, (br_l_seq r+1) `div` 2, br_n_cigar_op r * 4 ] 
  where
    r = BamRaw o s 
    good | S.length s < 32 = False
         | otherwise       = S.length s >= m
    m = sum [ 32, br_l_read_name r, br_l_seq r, (br_l_seq r+1) `div` 2, br_n_cigar_op r * 4 ] 

-- | Accessor for raw bam.
br_qname :: BamRaw -> S.ByteString
br_qname r@(BamRaw _ raw) = S.unsafeTake (br_l_read_name r) $ S.unsafeDrop 32 raw

br_l_read_name :: BamRaw -> Int  
br_l_read_name (BamRaw _ raw) = fromIntegral $ S.unsafeIndex raw 8 - 1

br_l_seq :: BamRaw -> Int
br_l_seq (BamRaw _ raw) = getInt raw 16

getInt16 :: (Num a, Bits a) => S.ByteString -> Int -> a
getInt16 s o = fromIntegral (S.unsafeIndex s o) .|. fromIntegral (S.unsafeIndex s $ o+1) `shiftL`  8

getInt :: (Num a, Bits a) => S.ByteString -> Int -> a
getInt s o = fromIntegral (S.unsafeIndex s $ o+0)             .|. fromIntegral (S.unsafeIndex s $ o+1) `shiftL`  8 .|.
             fromIntegral (S.unsafeIndex s $ o+2) `shiftL` 16 .|. fromIntegral (S.unsafeIndex s $ o+3) `shiftL` 24

br_n_cigar_op :: BamRaw -> Int
br_n_cigar_op (BamRaw _ raw) = getInt16 raw 12

br_flag :: BamRaw -> Int
br_flag (BamRaw _ raw) = getInt16 raw 14

br_isPaired, br_isProperlyPaired, br_isUnmapped, br_isMateUnmapped, br_isReversed,
    br_isMateReversed, br_isFirstMate, br_isSecondMate, br_isAuxillary, br_isFailsQC,
    br_isDuplicate :: BamRaw -> Bool

br_isPaired         = flip testBit  0 . br_flag
br_isProperlyPaired = flip testBit  1 . br_flag
br_isUnmapped       = flip testBit  2 . br_flag
br_isMateUnmapped   = flip testBit  3 . br_flag
br_isReversed       = flip testBit  4 . br_flag
br_isMateReversed   = flip testBit  5 . br_flag
br_isFirstMate      = flip testBit  6 . br_flag
br_isSecondMate     = flip testBit  7 . br_flag
br_isAuxillary      = flip testBit  8 . br_flag
br_isFailsQC        = flip testBit  9 . br_flag
br_isDuplicate      = flip testBit 10 . br_flag

br_rname :: BamRaw -> Refseq
br_rname (BamRaw _ raw) = Refseq $ getInt raw 0

br_pos :: BamRaw -> Int
br_pos (BamRaw _ raw) = getInt raw 4

br_mrnm :: BamRaw -> Refseq
br_mrnm (BamRaw _ raw) = Refseq $ getInt raw 20

br_mpos :: BamRaw -> Int
br_mpos (BamRaw _ raw) = getInt raw 24

br_isize :: BamRaw -> Int
br_isize (BamRaw _ raw) | i >= 0x80000000 = i - 0x100000000
                        | otherwise       = i
    where i :: Int
          i = getInt raw 28

br_aln_length :: BamRaw -> Int
br_aln_length br@(BamRaw _ raw) 
    | ncig == 0 = br_l_seq br
    | otherwise = sum [ x `shiftR` 4 | x <- map (getInt raw) $ take ncig [cig_off, cig_off+4 ..]
                                     , x .&. 0xF == 0 || x .&. 0xF == 2 || x .&. 0xF == 3 ]
    where
        !ncig    = br_n_cigar_op br
        !cig_off = 33 + br_l_read_name br

-- | Decode a BAM stream into raw entries.  Note that the entries can be
-- unpacked using @decodeBamEntry@.  Also note that this is an
-- Enumeratee in spirit, only the @BamMeta@ and @Refs@ need to get
-- passed separately.
decodeBam :: Monad m => (BamMeta -> Iteratee [BamRaw] m a) -> Iteratee Block m (Iteratee [BamRaw] m a)
decodeBam inner = do meta <- liftBlock get_bam_header
                     refs <- liftBlock get_ref_array
                     decodeBamLoop $ inner $! merge meta refs
  where
    get_bam_header  = do magic <- heads "BAM\SOH"
                         when (magic /= 4) $ fail "BAM signature not found"
                         hdr_len <- endianRead4 LSB
                         joinI $ i'take (fromIntegral hdr_len) $ parserToIteratee parseBamMeta

    get_ref_array = do nref <- endianRead4 LSB
                       foldM (\acc _ -> do 
                           nm <- endianRead4 LSB >>= i'getString . fromIntegral
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


decodeBamLoop :: Monad m => Enumeratee Block [BamRaw] m a
decodeBamLoop = eneeCheckIfDone loop
  where
    loop k = isFinished >>= loop' k
    loop' k True = return $ liftI k
    loop' k False = do off <- getOffset
                       raw <- liftBlock $ do
                                bsize <- endianRead4 LSB
                                when (bsize < 32) $ fail "short BAM record"
                                i'getString (fromIntegral bsize)
                       eneeCheckIfDone loop . k $ Chunk [bamRaw off raw]

-- | Stop gap solution for a cheap index.  We only get the first offset
-- from the linear index, which allows us to navigate to a target
-- sequence.  Will do the rest when I need it.
type BamIndex = UArray Refseq Int64

readBamIndex :: FilePath -> IO BamIndex
readBamIndex = fileDriver readBamIndex'

readBamIndex' :: MonadIO m => Iteratee S.ByteString m BamIndex
readBamIndex' = do magic <- heads "BAI\1"
                   when (magic /= 4) $ fail "BAI signature not found"
                   nref <- fromIntegral `liftM` endianRead4 LSB
                   if nref < 1 then return (array (toEnum 1, toEnum 0) [])
                               else get_array nref
  where
    get_array nref = do 
        arr <- liftIO $ ( newArray_ (toEnum 0, toEnum (nref-1)) :: IO (IOUArray Refseq Int64) )
        forM_ [toEnum 0 .. toEnum (nref-1)] $ \r -> do
            nbins <- fromIntegral `liftM` endianRead4 LSB
            replicateM_ nbins $ do
                _bin <- endianRead4 LSB -- "distinct bin", whatever that means
                nchunks <- fromIntegral `liftM` endianRead4 LSB
                replicateM_ nchunks $ endianRead8 LSB >> endianRead8 LSB

            nintv <- endianRead4 LSB
            o <- let loop acc 0 = return acc
                     loop acc n = do oo <- fromIntegral `liftM` endianRead8 LSB
                                     let !acc' = if oo == 0 then acc else min acc oo
                                     loop acc' (n-1)
                 in loop maxBound nintv                    
            liftIO $ writeArray arr r o 
        liftIO $ unsafeFreeze arr


writeRawBamFile :: MonadCatchIO m => FilePath -> BamMeta -> Iteratee [BamRaw] m ()
writeRawBamFile fp meta = 
    CIO.bracket (liftIO $ openBinaryFile fp WriteMode)
                (liftIO . hClose)
                (flip writeRawBamHandle meta) 

pipeRawBamOutput :: MonadIO m => BamMeta -> Iteratee [BamRaw] m ()
pipeRawBamOutput meta = encodeBamUncompressed meta =$ mapChunksM_ (liftIO . S.hPut stdout)

writeRawBamHandle :: MonadIO m => Handle -> BamMeta -> Iteratee [BamRaw] m ()
writeRawBamHandle hdl meta = encodeBam meta =$ mapChunksM_ (liftIO . S.hPut hdl)


mutateBamRaw :: BamRaw -> Mutator () -> BamRaw
mutateBamRaw (BamRaw vo br) mut = unsafePerformIO $ do
        S.useAsCStringLen br $ \(p,l) -> do
            runMutator mut p l
            bamRaw vo `fmap` S.packCStringLen (p,l)

newtype Mutator a = Mut { runMutator :: CString -> Int -> IO a }

instance Monad Mutator where
    return a = Mut $ \_ _ -> return a
    m >>= k  = Mut $ \p l -> runMutator m p l >>= \a -> runMutator (k a) p l

setFlag, setMpos, setIsize :: Int -> Mutator ()
setFlag  f = Mut $ \p _ -> pokeInt16 p 14 f
setMpos  x = Mut $ \p _ -> pokeInt32 p 24 x
setIsize x = Mut $ \p _ -> pokeInt32 p 28 x

setMrnm :: Refseq -> Mutator ()
setMrnm r = Mut $ \p _ -> pokeInt32 p 20 (unRefseq r)

pokeInt16 :: (Bits a, Integral a) => CString -> Int -> a -> IO ()
pokeInt16 p o x = do pokeElemOff p  o    . fromIntegral $        x   .&. 0xff
                     pokeElemOff p (o+1) . fromIntegral $ shiftR x 8 .&. 0xff

pokeInt32 :: (Bits a, Integral a) => CString -> Int -> a -> IO ()
pokeInt32 p o x = do pokeElemOff p  o    . fromIntegral $        x    .&. 0xff
                     pokeElemOff p (o+1) . fromIntegral $ shiftR x  8 .&. 0xff
                     pokeElemOff p (o+2) . fromIntegral $ shiftR x 16 .&. 0xff
                     pokeElemOff p (o+3) . fromIntegral $ shiftR x 24 .&. 0xff

