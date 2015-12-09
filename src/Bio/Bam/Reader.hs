{-# LANGUAGE BangPatterns, OverloadedStrings, FlexibleContexts #-}
module Bio.Bam.Reader (
    Block(..),
    decompressBgzfBlocks,
    decompressBgzf,
    compressBgzf,

    decodeBam,
    getBamRaw,
    decodeAnyBam,
    decodeAnyBamFile,

    BamrawEnumeratee,
    BamEnumeratee,
    isBamOrSam,

    isBam,
    isPlainBam,
    isGzipBam,
    isBgzfBam,

    decodeSam,
    decodeSam',

    decodeAnyBamOrSam,
    decodeAnyBamOrSamFile,

    concatInputs,
    concatDefaultInputs,
    mergeInputs,
    mergeDefaultInputs,
    combineCoordinates,
    combineNames,
                      ) where

import Bio.Base
import Bio.Bam.Header
import Bio.Bam.Rec
import Bio.Iteratee
import Bio.Iteratee.Bgzf
import Bio.Iteratee.ZLib hiding ( CompressionLevel )

import Control.Applicative
import Control.Arrow                ( (&&&) )
import Control.Monad
import Data.Attoparsec.ByteString   ( anyWord8 )
import Data.Char                    ( digitToInt )
import Data.Sequence                ( (|>) )
import Data.String                  ( fromString )
import System.Environment           ( getArgs )

import qualified Data.Attoparsec.ByteString.Char8   as P
import qualified Data.ByteString                    as B
import qualified Data.ByteString.Char8              as S
import qualified Data.Foldable                      as F
import qualified Data.Map.Strict                    as M
import qualified Data.Sequence                      as Z
import qualified Data.Vector.Generic                as V
import qualified Data.Vector.Storable               as VS
import qualified Data.Vector.Unboxed                as U

-- ^ Parsers for BAM and SAM.  We employ an @Iteratee@ interface, and we
-- strive to support everything possible in BAM.  The implementation of
-- nucleotides is somewhat lacking:  the "=" symbol is not understood.
--
-- TONOTDO:
-- - Reader for gzipped/bzipped/bgzf'ed SAM.  Storing SAM is a bad idea,
--   so why would anyone ever want to compress, much less index it?

type ByteString = B.ByteString
type BamrawEnumeratee m b = Enumeratee' BamMeta S.ByteString [BamRaw] m b
type BamEnumeratee m b = Enumeratee' BamMeta ByteString [BamRec] m b

isBamOrSam :: MonadIO m => Iteratee ByteString m (BamEnumeratee m a)
isBamOrSam = maybe decodeSam wrap `liftM` isBam
  where
    wrap enee it' = enee (\hdr -> mapStream unpackBam (it' hdr)) >>= lift . run


-- | Checks if a file contains BAM in any of the common forms,
-- then decompresses it appropriately.  If the stream doesn't contain
-- BAM at all, it is instead decoded as SAM.  Since SAM is next to
-- impossible to recognize reliably, we don't even try.  Any old junk is
-- decoded as SAM and will fail later.
decodeAnyBamOrSam :: MonadIO m => BamEnumeratee m a
decodeAnyBamOrSam it = isBamOrSam >>= \k -> k it

decodeAnyBamOrSamFile :: (MonadIO m, MonadMask m)
                      => FilePath -> (BamMeta -> Iteratee [BamRec] m a) -> m (Iteratee [BamRec] m a)
decodeAnyBamOrSamFile fn k = enumFileRandom defaultBufSize fn (decodeAnyBamOrSam k) >>= run

-- | Iteratee-style parser for SAM files, designed to be compatible with
-- the BAM parsers.  Parses plain uncompressed SAM, nothing else.  Since
-- it is supposed to work the same way as the BAM parser, it requires
-- the presense of the SQ header lines.  These are stripped from the
-- header text and turned into the symbol table.
decodeSam :: Monad m => (BamMeta -> Iteratee [BamRec] m a) -> Iteratee ByteString m (Iteratee [BamRec] m a)
decodeSam inner = joinI $ enumLinesBS $ do
    let pHeaderLine acc str = case P.parseOnly parseBamMetaLine str of Right f -> return $ f : acc
                                                                       Left e  -> fail $ e ++ ", " ++ show str
    meta <- liftM (foldr ($) mempty . reverse) (joinI $ breakE (not . S.isPrefixOf "@") $ foldStreamM pHeaderLine [])
    decodeSamLoop (meta_refs meta) (inner meta)

decodeSamLoop :: Monad m => Refs -> Enumeratee [ByteString] [BamRec] m a
decodeSamLoop refs inner = convStream (liftI parse_record) inner
  where !refs' = M.fromList $ zip [ nm | BamSQ { sq_name = nm } <- F.toList refs ] [toEnum 0..]
        ref x = M.findWithDefault invalidRefseq x refs'

        parse_record (EOF x) = icont parse_record x
        parse_record (Chunk []) = liftI parse_record
        parse_record (Chunk (l:ls)) | "@" `S.isPrefixOf` l = parse_record (Chunk ls)
        parse_record (Chunk (l:ls)) = case P.parseOnly (parseSamRec ref) l of
            Right  r -> idone [r] (Chunk ls)
            Left err -> icont parse_record (Just $ iterStrExc $ err ++ ", " ++ show l)

-- | Parser for SAM that doesn't look for a header.  Has the advantage
-- that it doesn't stall on a pipe that never delivers data.  Has the
-- disadvantage that it never reads the header and therefore needs a
-- list of allowed RNAMEs.
decodeSam' :: Monad m => Refs -> Enumeratee ByteString [BamRec] m a
decodeSam' refs inner = joinI $ enumLinesBS $ decodeSamLoop refs inner

parseSamRec :: (ByteString -> Refseq) -> P.Parser BamRec
parseSamRec ref = mkBamRec
                  <$> word <*> num <*> (ref <$> word) <*> (subtract 1 <$> num)
                  <*> (Q <$> num') <*> (VS.fromList <$> cigar) <*> rnext <*> (subtract 1 <$> num)
                  <*> snum <*> sequ <*> quals <*> exts <*> pure 0
  where
    sep      = P.endOfInput <|> () <$ P.char '\t'
    word     = P.takeTill ((==) '\t') <* sep
    num      = P.decimal <* sep
    num'     = P.decimal <* sep
    snum     = P.signed P.decimal <* sep

    rnext    = id <$ P.char '=' <* sep <|> const . ref <$> word
    sequ     = {-# SCC "parseSamRec/sequ" #-}
               (V.empty <$ P.char '*' <|>
               V.fromList . map toNucleotides . S.unpack <$> P.takeWhile is_nuc) <* sep

    quals    = {-# SCC "parseSamRec/quals" #-} defaultQs <$ P.char '*' <* sep <|> bsToVec <$> word
        where
            defaultQs sq = VS.replicate (V.length sq) (Q 0xff)
            bsToVec qs _ = VS.fromList . map (Q . subtract 33) $ B.unpack qs

    cigar    = [] <$ P.char '*' <* sep <|>
               P.manyTill (flip (:*) <$> P.decimal <*> cigop) sep

    cigop    = P.choice $ zipWith (\c r -> r <$ P.char c) "MIDNSHP" [Mat,Ins,Del,Nop,SMa,HMa,Pad]
    exts     = ext `P.sepBy` sep
    ext      = (\a b v -> (fromString [a,b],v)) <$> P.anyChar <*> P.anyChar <*> (P.char ':' *> value)

    value    = P.char 'A' *> P.char ':' *> (Char <$>               anyWord8) <|>
               P.char 'i' *> P.char ':' *> (Int  <$>     P.signed P.decimal) <|>
               P.char 'Z' *> P.char ':' *> (Text <$> P.takeTill ((==) '\t')) <|>
               P.char 'H' *> P.char ':' *> (Bin  <$>               hexarray) <|>
               P.char 'f' *> P.char ':' *> (Float . realToFrac <$> P.double) <|>
               P.char 'B' *> P.char ':' *> (
                    P.satisfy (P.inClass "cCsSiI") *> (intArr   <$> many (P.char ',' *> P.signed P.decimal)) <|>
                    P.char 'f'                     *> (floatArr <$> many (P.char ',' *> P.double)))

    intArr   is = IntArr   $ U.fromList is
    floatArr fs = FloatArr $ U.fromList $ map realToFrac fs
    hexarray    = B.pack . repack . S.unpack <$> P.takeWhile (P.inClass "0-9A-Fa-f")
    repack (a:b:cs) = fromIntegral (digitToInt a * 16 + digitToInt b) : repack cs ; repack _ = []
    is_nuc = P.inClass "acgtswkmrybdhvnACGTSWKMRYBDHVN"

    mkBamRec nm fl rn po mq cg rn' mp is sq qs' =
                BamRec nm fl rn po mq cg (rn' rn) mp is sq (qs' sq)

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
  where u ? v = if (b_rname &&& b_pos) (unpackBam u) < (b_rname &&& b_pos) (unpackBam v) then Less else NotLess

{-# INLINE combineNames #-}
combineNames :: Monad m => BamMeta -> Enumeratee [BamRaw] [BamRaw] (Iteratee [BamRaw] m) a
combineNames _ = mergeSortStreams (?)
  where u ? v = case b_qname (unpackBam u) `compareNames` b_qname (unpackBam v) of LT -> Less ; _ -> NotLess

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
