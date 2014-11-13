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

    CigOp(..),
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
    br_extflag,

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
    br_bin,
    br_mapq,
    br_mrnm,
    br_mpos,
    br_isize,
    br_seq_at,
    br_qual_at,
    br_cigar_at,

    br_findExtension,
    br_extAsInt,
    br_extAsString,

    Mutator,
    mutateBamRaw,
    removeExt,
    appendStringExt,
    setPos,
    setBin,
    setMapq,
    setFlag,
    setMrnm,
    setMpos,
    setIsize
) where

import Bio.Base
import Bio.Bam.Header
import Bio.Iteratee
import Bio.Iteratee.ZLib hiding ( CompressionLevel )
import Bio.Iteratee.Bgzf

import Control.Monad
import Data.Binary.Put
import Data.Bits                    ( Bits, shiftL, shiftR, (.&.), (.|.), testBit )
import Data.Int                     ( Int32, Int16, Int8 )
import Data.Monoid
import Data.Sequence                ( (|>) )
import Data.Word                    ( Word32, Word16 )
import Foreign.C.String             ( CString )
import Foreign.Ptr                  ( plusPtr )
import Foreign.Marshal.Utils        ( moveBytes )
import Foreign.Storable             ( pokeElemOff )
import System.Environment           ( getArgs )
import System.IO
import System.IO.Unsafe

import qualified Control.Monad.Catch            as C
import qualified Data.ByteString                as S
import qualified Data.ByteString.Char8          as SC
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
-- recognized as Bam, a decoder (suitable Enumeratee) for it is returned.
isBam, isEmptyBam, isPlainBam, isBgzfBam, isGzipBam :: MonadIO m
    => Iteratee S.ByteString m (Maybe (BamrawEnumeratee m a))
isBam = firstOf [ isEmptyBam, isPlainBam, isBgzfBam, isGzipBam ]
  where
    firstOf [] = return Nothing
    firstOf (k:ks) = k >>= maybe (firstOf ks) (return . Just)

isEmptyBam = (\e -> if e then Just (\k -> return $ k mempty) else Nothing) `liftM` isFinished

isPlainBam = (\n -> if n == 4 then Just (joinI . decompressPlain . decodeBam) else Nothing)
             `liftM` iLookAhead (heads "BAM\SOH")

-- Interesting... iLookAhead interacts badly with the parallel
-- decompression of BGZF.  (The chosen interface doesn't allow the EOF
-- signal to be passed on.)  One workaround would be to run sequential
-- BGZF decompression to check if the content is BAM, but since BGZF is
-- actually GZip in disguise, the easier workaround if to use the
-- ordinary GZip decompressor.
-- (A clean workaround would be an @Alternative@ instance for
-- @Iteratee@.)
isBgzfBam  = do b <- isBgzf
                k <- if b then iLookAhead $ joinI $ enumInflate GZip defaultDecompressParams isPlainBam
                          else return Nothing
                return $ (\_ -> (joinI . decompressBgzfBlocks . decodeBam)) `fmap` k

isGzipBam  = do b <- isGzip
                k <- if b then iLookAhead $ joinI $ enumInflate GZip defaultDecompressParams isPlainBam
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
    eneeBam  = eneeCheckIfDone (\k -> eneeBam2 . k $ Chunk header)
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
{-# INLINE br_qname #-}
br_qname :: BamRaw -> S.ByteString
br_qname r@(BamRaw _ raw) = S.unsafeTake (br_l_read_name r) $ S.unsafeDrop 32 raw

{-# INLINE br_l_read_name #-}
br_l_read_name :: BamRaw -> Int
br_l_read_name (BamRaw _ raw) = fromIntegral $ S.unsafeIndex raw 8 - 1

{-# INLINE br_l_seq #-}
br_l_seq :: BamRaw -> Int
br_l_seq (BamRaw _ raw) = getInt raw 16

{-# INLINE getInt16 #-}
getInt16 :: (Num a, Bits a) => S.ByteString -> Int -> a
getInt16 s o = fromIntegral (S.unsafeIndex s o) .|. fromIntegral (S.unsafeIndex s $ o+1) `shiftL`  8

-- | Load an unaligned, little-endian int.  This is probably quite slow
-- and unnecessary on some platforms.  On i386, ix86_64 and powerPC, we
-- could cast the pointer and do a direct load.  Other may have special
-- primitives.  Worth investigating?
{-# INLINE getInt #-}
getInt :: (Num a, Bits a) => S.ByteString -> Int -> a
getInt s o = fromIntegral (S.unsafeIndex s $ o+0)             .|. fromIntegral (S.unsafeIndex s $ o+1) `shiftL`  8 .|.
             fromIntegral (S.unsafeIndex s $ o+2) `shiftL` 16 .|. fromIntegral (S.unsafeIndex s $ o+3) `shiftL` 24

{-# INLINE br_n_cigar_op #-}
br_n_cigar_op :: BamRaw -> Int
br_n_cigar_op (BamRaw _ raw) = getInt16 raw 12

{-# INLINE br_flag #-}
br_flag :: BamRaw -> Int
br_flag (BamRaw _ raw) = getInt16 raw 14

{-# INLINE br_extflag #-}
br_extflag :: BamRaw -> Int
br_extflag br = shiftL ef 16 .|. ff
  where !ff = br_flag br
        !ef = br_extAsInt (br_extAsInt 0 "XF" br) "FF" br

{-# INLINE br_isPaired         #-}
{-# INLINE br_isProperlyPaired #-}
{-# INLINE br_isUnmapped       #-}
{-# INLINE br_isMateUnmapped   #-}
{-# INLINE br_isReversed       #-}
{-# INLINE br_isMateReversed   #-}
{-# INLINE br_isFirstMate      #-}
{-# INLINE br_isSecondMate     #-}
{-# INLINE br_isAuxillary      #-}
{-# INLINE br_isFailsQC        #-}
{-# INLINE br_isDuplicate      #-}
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

{-# INLINE br_rname #-}
br_rname :: BamRaw -> Refseq
br_rname (BamRaw _ raw) = Refseq $ getInt raw 0

{-# INLINE br_bin #-}
br_bin :: BamRaw -> Int
br_bin (BamRaw _ raw) = getInt16 raw 10

{-# INLINE br_mapq #-}
br_mapq :: BamRaw -> Qual
br_mapq (BamRaw _ raw) = Q $ S.unsafeIndex raw 9

{-# INLINE br_pos #-}
br_pos :: BamRaw -> Int
br_pos (BamRaw _ raw) = getInt raw 4

{-# INLINE br_mrnm #-}
br_mrnm :: BamRaw -> Refseq
br_mrnm (BamRaw _ raw) = Refseq $ getInt raw 20

{-# INLINE br_mpos #-}
br_mpos :: BamRaw -> Int
br_mpos (BamRaw _ raw) = getInt raw 24

{-# INLINE br_isize #-}
br_isize :: BamRaw -> Int
br_isize (BamRaw _ raw) | i >= 0x80000000 = i - 0x100000000
                        | otherwise       = i
    where i :: Int
          i = getInt raw 28

{-# INLINE br_seq_at #-}
br_seq_at :: BamRaw -> Int -> Nucleotides
br_seq_at br@(BamRaw _ raw) i
    | even    i = Ns $ (S.unsafeIndex raw (off0 + i `div` 2) `shiftR` 4) .&. 0xF
    | otherwise = Ns $  S.unsafeIndex raw (off0 + i `div` 2)             .&. 0xF
  where
    off0 = sum [ 33, br_l_read_name br, 4 * br_n_cigar_op br ]

{-# INLINE br_qual_at #-}
br_qual_at :: BamRaw -> Int -> Qual
br_qual_at br@(BamRaw _ raw) i = Q $ S.unsafeIndex raw (off0 + i)
  where
    off0 = sum [ 33, br_l_read_name br, 4 * br_n_cigar_op br, (br_l_seq br + 1) `div` 2]

{-# INLINE br_cigar_at #-}
br_cigar_at :: BamRaw -> Int -> (CigOp, Int)
br_cigar_at br@(BamRaw _ raw) i = (co,cl)
  where
    !cig_off = 33 + br_l_read_name br
    !cig_val = getInt raw $ cig_off + i * 4
    !cl = fromIntegral cig_val `shiftR` 4
    !co = case cig_val .&. 0xf :: Int of {
                1 -> Ins ; 2 -> Del ;
                3 -> Nop ; 4 -> SMa ; 5 -> HMa ;
                6 -> Pad ; _ -> Mat }

{-# INLINE br_aln_length #-}
br_aln_length :: BamRaw -> Int
br_aln_length br@(BamRaw _ raw)
    = sum [ x `shiftR` 4 | x <- map (getInt raw) $ take ncig [cig_off, cig_off+4 ..]
          , x .&. 0xF == 0 || x .&. 0xF == 2 || x .&. 0xF == 3 ]
  where
    !ncig    = br_n_cigar_op br
    !cig_off = 33 + br_l_read_name br

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


mutateBamRaw :: BamRaw -> Mutator () -> BamRaw
mutateBamRaw (BamRaw vo br) mut = unsafePerformIO $ do
        S.useAsCStringLen br $ \(p,l) -> do
            (l',frags, ()) <- runMutator mut p l []
            if l' <= l then do f1 <- S.packCStringLen (p,l')
                               return $! bamRaw vo $! S.concat (f1 : reverse frags)
                       else error "broken Mutator: length must never increase"

newtype Mutator a = Mut { runMutator :: CString -> Int -> [S.ByteString] -> IO (Int,[S.ByteString], a) }

instance Monad Mutator where
    {-# INLINE return #-}
    return a = Mut $ \_ l fs -> return (l,fs,a)
    {-# INLINE (>>=) #-}
    m >>= k  = Mut $ \p l fs -> runMutator m p l fs >>= \(l',fs',a) -> runMutator (k a) p l' fs'

{-# INLINE passL #-}
passL :: IO a -> Int -> [S.ByteString] -> IO (Int,[S.ByteString],a)
passL io = \l fs -> io >>= \a -> return (l,fs,a)

{-# INLINE setFlag  #-}
{-# INLINE setPos  #-}
{-# INLINE setMpos  #-}
{-# INLINE setIsize #-}
setPos, setFlag, setMpos, setIsize :: Int -> Mutator ()
setPos   x = Mut $ \p -> passL $ pokeInt32 p  4 x
setFlag  f = Mut $ \p -> passL $ pokeInt16 p 14 f
setMpos  x = Mut $ \p -> passL $ pokeInt32 p 24 x
setIsize x = Mut $ \p -> passL $ pokeInt32 p 28 x

{-# INLINE setMapq #-}
setMapq :: Word8 -> Mutator ()
setMapq q = Mut $ \p -> passL $ pokeInt8 p 9 q

{-# INLINE setBin #-}
setBin :: Int -> Int -> Mutator ()
setBin b l = Mut $ \p -> passL $ pokeInt16 p 10 (distinctBin b l)

{-# INLINE setMrnm #-}
setMrnm :: Refseq -> Mutator ()
setMrnm r = Mut $ \p -> passL $ pokeInt32 p 20 (unRefseq r)

{-# INLINE pokeInt8 #-}
pokeInt8 :: Integral a => CString -> Int -> a -> IO ()
pokeInt8 p o = pokeElemOff p o . fromIntegral

{-# INLINE pokeInt16 #-}
pokeInt16 :: (Bits a, Integral a) => CString -> Int -> a -> IO ()
pokeInt16 p o x = do pokeElemOff p  o    . fromIntegral $        x   .&. 0xff
                     pokeElemOff p (o+1) . fromIntegral $ shiftR x 8 .&. 0xff

{-# INLINE pokeInt32 #-}
pokeInt32 :: (Bits a, Integral a) => CString -> Int -> a -> IO ()
pokeInt32 p o x = do pokeElemOff p  o    . fromIntegral $        x    .&. 0xff
                     pokeElemOff p (o+1) . fromIntegral $ shiftR x  8 .&. 0xff
                     pokeElemOff p (o+2) . fromIntegral $ shiftR x 16 .&. 0xff
                     pokeElemOff p (o+3) . fromIntegral $ shiftR x 24 .&. 0xff



-- Find an extension field, return offset in BamRaw data.
{-# INLINE br_findExtension #-}
br_findExtension :: String -> BamRaw -> Maybe (Int,Int,Int)
br_findExtension [u,v] br@(BamRaw _ r) = go off0
  where
    off0 = sum [ 33, br_l_read_name br, 4 * br_n_cigar_op br, br_l_seq br, (br_l_seq br +1) `div` 2 ]
    go !o | o >= S.length r - 3                        = Nothing
          | SC.index r o == u && SC.index r (o+1) == v = Just (o, o+2, skip o)
          | otherwise                                  = go (skip o)

    skip !o = case SC.index r (o+2) of
        'Z' -> skipNul $ o + 3
        'H' -> skipNul $ o + 3
        'B' -> o + 7 + sizeof (SC.index r (o+3)) * getInt r (o+4)
        xxx -> o + 3 + sizeof xxx

    skipNul !o | S.length r  == o = o
               | S.index r o == 0 = o+1
               | otherwise        = skipNul (o+1)

    sizeof 'A' = 1
    sizeof 'c' = 1
    sizeof 'C' = 1
    sizeof 's' = 2
    sizeof 'S' = 2
    sizeof 'i' = 4
    sizeof 'I' = 4
    sizeof 'f' = 4
    sizeof  x  = error $ "unknown fields type: " ++ show x
br_findExtension _ _ = error "illegal key, must be two characters"

{-# INLINE br_extAsInt #-}
br_extAsInt :: Int -> String -> BamRaw -> Int
br_extAsInt d k br@(BamRaw _ r) = case br_findExtension k br of
        Just (_,o,_) | SC.index r o == 'c' -> fromIntegral               (S.index r (o+1))
                     | SC.index r o == 'C' -> fromIntegral (fromIntegral (S.index r (o+1)) :: Int8)
                     | SC.index r o == 's' -> fromIntegral (getInt16 r (o+1) :: Int16)
                     | SC.index r o == 'S' -> fromIntegral (getInt16 r (o+1) :: Word16)
                     | SC.index r o == 'i' -> fromIntegral (getInt   r (o+1) :: Int32)
                     | SC.index r o == 'I' -> fromIntegral (getInt   r (o+1) :: Word32)
        _                                  -> d

{-# INLINE br_extAsString #-}
br_extAsString :: String -> BamRaw -> S.ByteString
br_extAsString k br@(BamRaw _ r) = case br_findExtension k br of
        Just (_,o,_) | SC.index r o == 'A' -> S.singleton (S.index r (o+1))
                     | SC.index r o == 'Z' -> S.takeWhile (/= 0) $ S.drop (o+1) r
                     | SC.index r o == 'H' -> S.takeWhile (/= 0) $ S.drop (o+1) r
        _                                  -> S.empty

{-# INLINE removeExt #-}
removeExt :: String -> Mutator ()
removeExt key = Mut $ \p l fs -> do
    r <- S.unsafePackCStringLen (p,l)
    case br_findExtension key (bamRaw 0 r) of
        Nothing      -> return (l,fs,())
        Just (a,_,b) -> do moveBytes (p `plusPtr` a) (p `plusPtr` b) (l-b)
                           return $ (l-(b-a), fs, ())

{-# INLINE appendStringExt #-}
appendStringExt :: String -> S.ByteString -> Mutator ()
appendStringExt [u,v] value = Mut $ \_ l fs -> return (l,f:fs,())
  where
    f = S.concat [ SC.singleton u, SC.singleton v, SC.singleton 'Z', value, S.singleton 0 ]
appendStringExt _ _ = error "illegal key, must be two characters"

