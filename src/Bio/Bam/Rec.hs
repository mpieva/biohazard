{-# LANGUAGE OverloadedStrings, PatternGuards, BangPatterns, NoMonomorphismRestriction #-}

-- TODO:
-- - Automatic creation of some kind of index.  If possible, this should
--   be the standard index for sorted BAM and/or the newer CSI format.
--   Optionally, a block index for slicing of large files, even unsorted
--   ones.  Maybe an index by name and an index for group-sorted files.
--   Sensible indices should be generated whenever a file is written.
-- - Same for statistics.  Something like "flagstats" could always be
--   written.  Actually, having @writeBamHandle@ return enhanced
--   flagstats as a result might be even better.
--
-- TONOTDO:
-- - Reader for gzipped/bzipped/bgzf'ed SAM.  Storing SAM is a bad idea,
--   so why would anyone ever want to compress, much less index it?

module Bio.Bam.Rec (
    Block,
    BamEnumeratee,
    isBamOrSam,

    decodeBamEntry,
    encodeBamEntry,
    encodeSamEntry,

    decodeSam,
    decodeSam',

    decodeAnyBamOrSam,
    decodeAnyBamOrSamFile,

    writeBamFile,
    writeBamHandle,
    pipeBamOutput,
    pipeRawSamOutput,
    pipeSamOutput,

    BamRec(..),
    nullBamRec,
    getMd,

    Nucleotides(..),
    Extensions, Ext(..),
    extAsInt, extAsString, setQualFlag,

    isPaired,
    isProperlyPaired,
    isUnmapped,
    isMateUnmapped,
    isReversed,
    isMateReversed,
    isFirstMate,
    isSecondMate,
    isAuxillary,
    isFailsQC,
    isDuplicate,
    isTrimmed,
    isMerged,
    type_mask,

    Word32
) where

import Bio.Base
import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Iteratee

import Control.Monad
import Control.Applicative
import Data.Array.IArray
import Data.Array.Unboxed
import Data.Attoparsec.ByteString   ( anyWord8 )
import Data.Binary.Builder          ( toLazyByteString )
import Data.Binary.Get
import Data.Binary.Put
import Data.Bits                    ( testBit, shiftL, shiftR, (.&.), (.|.), complement )
import Data.ByteString              ( ByteString )
import Data.Char                    ( ord, digitToInt )
import Data.Int                     ( Int32 )
import Data.Monoid                  ( mempty )
import Data.Vector.Unboxed          ( (!?) )
import Data.Word                    ( Word32 )
import Foreign.Marshal.Alloc        ( alloca )
import Foreign.Ptr                  ( castPtr )
import Foreign.Storable             ( peek, poke )
import System.IO
import System.IO.Unsafe             ( unsafePerformIO )

import qualified Data.Attoparsec.ByteString.Char8   as P
import qualified Data.ByteString                    as B
import qualified Data.ByteString.Char8              as S
import qualified Data.ByteString.Lazy.Char8         as L
import qualified Data.Foldable                      as F
import qualified Data.Iteratee                      as I
import qualified Data.Map                           as M
import qualified Data.Vector.Unboxed                as V

-- ^ Parsers and Printers for BAM and SAM.  We employ an @Iteratee@
-- interface, and we strive to support everything possible in BAM.  So
-- far, the implementation of the nucleotides is somewhat lacking:  we
-- do not have support for ambiguity codes, and the "=" symbol is not
-- understood.

-- | internal representation of a BAM record
data BamRec = BamRec {
        b_qname :: {-# UNPACK #-} !Seqid,
        b_flag  :: {-# UNPACK #-} !Int,
        b_rname :: {-# UNPACK #-} !Refseq,
        b_pos   :: {-# UNPACK #-} !Int,
        b_mapq  :: {-# UNPACK #-} !Qual,
        b_cigar :: Cigar,
        b_mrnm  :: {-# UNPACK #-} !Refseq,
        b_mpos  :: {-# UNPACK #-} !Int,
        b_isize :: {-# UNPACK #-} !Int,
        b_seq   :: !(V.Vector Nucleotides),
        b_qual  :: !ByteString,         -- ^ quality, may be empty
        b_exts  :: Extensions,
        b_virtual_offset :: {-# UNPACK #-} !FileOffset -- ^ virtual offset for indexing purposes
    } deriving Show

nullBamRec :: BamRec
nullBamRec = BamRec {
        b_qname = S.empty,
        b_flag  = flagUnmapped,
        b_rname = invalidRefseq,
        b_pos   = invalidPos,
        b_mapq  = Q 0,
        b_cigar = Cigar [],
        b_mrnm  = invalidRefseq,
        b_mpos  = invalidPos,
        b_isize = 0,
        b_seq   = V.empty,
        b_qual  = S.empty,
        b_exts  = M.empty,
        b_virtual_offset = 0
    }

getMd :: BamRec -> Maybe [MdOp]
getMd r = case M.lookup "MD" $ b_exts r of
    Just (Text mdfield) -> readMd mdfield
    Just (Char mdfield) -> readMd $ B.singleton mdfield
    _                   -> Nothing

type BamEnumeratee m b = Enumeratee' BamMeta ByteString [BamRec] m b

isBamOrSam :: MonadIO m => Iteratee ByteString m (BamEnumeratee m a)
isBamOrSam = maybe decodeSam wrap `liftM` isBam
  where
    wrap enee it' = enee (\hdr -> I.mapStream decodeBamEntry (it' hdr)) >>= lift . run


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


-- | Decodes a raw block into a @BamRec@.
decodeBamEntry :: BamRaw -> BamRec
decodeBamEntry br = case pushEndOfInput $ runGetIncremental go `pushChunk` raw_data br of
        Fail _ _ m -> error m
        Partial  _ -> error "incomplete BAM record"
        Done _ _ r -> fixup_bam_rec r
  where
    go = do !rid       <- Refseq       <$> getWord32le
            !start     <- fromIntegral <$> getWord32le
            !namelen   <- fromIntegral <$> getWord8
            !mapq      <-            Q <$> getWord8
            !_bin      <-                  getWord16le
            !cigar_len <- fromIntegral <$> getWord16le
            !flag      <- fromIntegral <$> getWord16le
            !read_len  <- fromIntegral <$> getWord32le
            !mate_rid  <- Refseq       <$> getWord32le
            !mate_pos  <- fromIntegral <$> getWord32le
            !ins_size  <- fromIntegral <$> getWord32le
            !read_name <- S.init       <$> getByteString namelen
            !cigar     <- Cigar . map decodeCigar <$> replicateM cigar_len getWord32le
            !qry_seq   <- getByteString $ (read_len+1) `div` 2
            !qual <- (\qs -> if B.all (0xff ==) qs then B.empty else qs) <$> getByteString read_len
            !exts <- getExtensions M.empty

            return $ BamRec read_name flag rid start mapq cigar
                            mate_rid mate_pos ins_size
                            (V.fromListN read_len $ expand qry_seq)
                            qual exts (virt_offset br)

    expand t = if S.null t then [] else let x = B.head t in Ns (x `shiftR` 4) : Ns (x .&. 0xf) : expand (B.tail t)

    decodeCigar c | cc <= fromEnum (maxBound :: CigOp) = (toEnum cc, cl)
                  | otherwise = error "unknown Cigar operation"
      where cc = fromIntegral c .&. 0xf; cl = fromIntegral c `shiftR` 4

-- | fixes BAM records for changed conventions
fixup_bam_rec :: BamRec -> BamRec
fixup_bam_rec b =
    (if b_flag b .&. flagLowQuality /= 0 then setQualFlag 'Q' else id) $          -- low qual, new convention
    (if b_flag b .&. flagLowComplexity /= 0 then setQualFlag 'C' else id) $       -- low complexity, new convention
    b { b_flag = fixPP $ oflags .|. muflag .|. tflags .|. shiftL eflags 16        -- extended flags
      , b_exts = cleaned_exts }
  where
        -- removes old flag abuse
        flags' = b_flag b .&. complement (flagLowQuality .|. flagLowComplexity)
        oflags | flags' .&. flagPaired == 0 = flags' .&. complement (flagFirstMate .|. flagSecondMate)
               | otherwise                  = flags'

        -- set "mate unmapped" if self coordinates and mate coordinates are equal, but self is paired and mapped
        -- (BWA forgets this flag for invalid mate alignments)
        muflag = if mu then flagMateUnmapped else 0
        mu = and [ isPaired b, not (isUnmapped b)
                 , isReversed b == isMateReversed b
                 , b_rname b == b_mrnm b, b_pos b == b_mpos b ]

        -- merged & trimmed from old flag abuse
        is_merged  = flags' .&. (flagPaired .|. flagFirstMate .|. flagSecondMate) == flagFirstMate .|. flagSecondMate
        is_trimmed = flags' .&. (flagPaired .|. flagFirstMate .|. flagSecondMate) == flagSecondMate

        tflags = (if is_merged  then flagMerged  else 0) .|.
                 (if is_trimmed then flagTrimmed else 0)

        -- extended flags, renamed to avoid collision with BWA
        -- Goes like this:  if FF is there, use and remove it.  Else
        -- check if XF is there _and_is_numeric_.  If so, use it and
        -- remove it.  Else use 0 and leave it alone.  Note that this
        -- solves the collision with BWA, since BWA puts a character
        -- there, not an int.
        (eflags, cleaned_exts) = case (M.lookup "FF" (b_exts b), M.lookup "XF" (b_exts b)) of
                ( Just (Int i), _ ) -> (i, M.delete "FF" (b_exts b))
                ( _, Just (Int i) ) -> (i, M.delete "XF" (b_exts b))
                (       _,_       ) -> (0,                b_exts b )

        -- if either mate is unmapped, remove "properly paired"
        fixPP f | f .&. (flagUnmapped .|. flagMateUnmapped) == 0 = f
                | otherwise = f .&. complement flagProperlyPaired

        flagLowQuality    =  0x800
        flagLowComplexity = 0x1000

-- | A collection of extension fields.  The key is actually only two @Char@s, but that proved impractical.
-- (Hmm... we could introduce a Key type that is a 16 bit int, then give
-- it an @instance IsString@... practical?)
type Extensions = M.Map String Ext

data Ext = Int Int | Float Float | Text ByteString | Bin ByteString | Char Word8
         | IntArr (UArray Int Int) | FloatArr (UArray Int Float)
    deriving (Show, Eq, Ord)

getExtensions :: Extensions -> Get Extensions
getExtensions m = getExt <|> return m
  where
    getExt :: Get Extensions
    getExt = do
            key <- (\a b -> [w2c a, w2c b]) <$> getWord8 <*> getWord8
            typ <- getWord8
            let cont v = getExtensions $! M.insert key v m
            case w2c typ of
                    'Z' -> cont . Text =<< getByteStringNul
                    'H' -> cont . Bin  =<< getByteStringNul
                    'A' -> cont . Char =<< getWord8
                    'f' -> cont . Float . to_float =<< getWord32le
                    'B' -> do tp <- getWord8
                              n <- fromIntegral <$> getWord32le
                              case w2c tp of
                                 'f' -> cont . FloatArr . listArray (0,n) . map to_float =<< replicateM (n+1) getWord32le
                                 x | Just get <- M.lookup x get_some_int -> cont . IntArr . listArray (0,n) =<< replicateM (n+1) get
                                   | otherwise                           -> fail $ "array type code " ++ show x ++ " not recognized"
                    x | Just get <- M.lookup x get_some_int -> cont . Int =<< get
                      | otherwise                           -> fail $ "type code " ++ show x ++ " not recognized"

    to_float :: Word32 -> Float
    to_float word = unsafePerformIO $ alloca $ \buf ->
                    poke (castPtr buf) word >> peek buf

    get_some_int :: M.Map Char (Get Int)
    get_some_int = M.fromList $ zip "cCsSiI" [
                        fromIntegral <$> getWord8,
                        fromIntegral <$> getWord8,
                        fromIntegral <$> getWord16le,
                        fromIntegral <$> getWord16le,
                        fromIntegral <$> getWord32le,
                        fromIntegral <$> getWord32le ]



getByteStringNul :: Get ByteString
getByteStringNul = S.init <$> (lookAhead (get_len 1) >>= getByteString)
  where
    get_len l = getWord8 >>= \w -> if w == 0 then return l else get_len $! l+1




encodeBamEntry :: BamRec -> BamRaw
encodeBamEntry = bamRaw 0 . S.concat . L.toChunks . runPut . putEntry
  where
    putEntry  b = do putWord32le   $ unRefseq $ b_rname b
                     put_int_32    $ b_pos b
                     put_int_8     $ S.length (b_qname b) + 1
                     put_int_8     $ unQ (b_mapq b)
                     put_int_16    $ distinctBin (b_pos b) (cigarToAlnLen (b_cigar b))
                     put_int_16    $ length $ unCigar $ b_cigar b
                     put_int_16    $ b_flag b
                     put_int_32    $ V.length $ b_seq b
                     putWord32le   $ unRefseq $ b_mrnm b
                     put_int_32    $ b_mpos b
                     put_int_32    $ b_isize b
                     putByteString $ b_qname b
                     putWord8 0
                     mapM_ (put_int_32 . encodeCigar) $ unCigar $ b_cigar b
                     putSeq $ b_seq b
                     putByteString $ if not (S.null (b_qual b)) then b_qual b
                                     else B.replicate (V.length $ b_seq b) 0xff
                     forM_ (M.toList $ more_exts b) $ \(k,v) ->
                        case k of [c,d] -> putChr c >> putChr d >> putValue v
                                  _     -> error $ "invalid field key " ++ show k

    more_exts :: BamRec -> Extensions
    more_exts b = if xf /= 0 then x' else b_exts b
        where xf = b_flag b `shiftR` 16
              x' = M.insert "FF" (Int xf) $ b_exts b

    encodeCigar :: (CigOp,Int) -> Int
    encodeCigar (op,l) = fromEnum op .|. l `shiftL` 4

    putSeq :: V.Vector Nucleotides -> Put
    putSeq v = case v !? 0 of
                 Nothing -> return ()
                 Just a  -> case v !? 1 of
                    Nothing -> putWord8 (unNs a `shiftL` 4)
                    Just b  -> do putWord8 (unNs a `shiftL` 4 .|. unNs b)
                                  putSeq (V.drop 2 v)

-- | writes BAM encoded stuff to a @Handle@
-- We generate BAM with dynamic blocks, then stream them out to the file.
--
-- XXX This could write indexes on the side---a simple block index
-- for MapReduce style slicing, a standard BAM index or a name index
-- would be possible.
writeBamHandle :: Handle -> BamMeta -> Iteratee [BamRec] IO ()
writeBamHandle hdl meta = I.mapStream encodeBamEntry =$ writeRawBamHandle hdl meta

-- | writes BAM encoded stuff to a file
-- XXX This should(!) write indexes on the side---a simple block index
-- for MapReduce style slicing, a standard BAM index or a name index
-- would be possible.  When writing to a file, this makes even more
-- sense than when writing to a @Handle@.
writeBamFile :: FilePath -> BamMeta -> Iteratee [BamRec] IO ()
writeBamFile fp meta = I.mapStream encodeBamEntry =$ writeRawBamFile fp meta

-- | write BAM encoded stuff to stdout
-- This send uncompressed BAM to stdout.  Useful for piping to other
-- tools.
pipeBamOutput :: BamMeta -> Iteratee [BamRec] IO ()
pipeBamOutput meta = I.mapStream encodeBamEntry =$ pipeRawBamOutput meta

-- | write in SAM format to stdout
-- This is useful for piping to other tools (say, AWK scripts) or for
-- debugging.  No convenience function to send SAM to a file exists,
-- because that's a stupid idea.
pipeSamOutput :: MonadIO m => BamMeta -> Iteratee [BamRec] m ()
pipeSamOutput meta = do liftIO . L.putStr . toLazyByteString $ showBamMeta meta
                        mapStreamM_ $ \b -> liftIO . putStr $ encodeSamEntry (meta_refs meta) b "\n"

pipeRawSamOutput :: MonadIO m => BamMeta -> Iteratee [BamRaw] m ()
pipeRawSamOutput hdr = joinI $ mapStream decodeBamEntry $ pipeSamOutput hdr

put_int_32, put_int_16, put_int_8 :: Integral a => a -> Put
put_int_32 = putWord32le . fromIntegral
put_int_16 = putWord16le . fromIntegral
put_int_8  = putWord8 . fromIntegral

putChr :: Char -> Put
putChr = putWord8 . fromIntegral . ord

putValue :: Ext -> Put
putValue v = case v of
    Text t      -> putChr 'Z' >> putByteString t >> putWord8 0
    Bin b       -> putChr 'H' >> putByteString b >> putWord8 0
    Char c      -> putChr 'A' >> putWord8 c
    Float f     -> putChr 'f' >> put_int_32 (fromFloat f)
    Int i       -> case put_some_int [i] of (c,op) -> putChr c >> op i
    FloatArr fa -> putChr 'B' >> putChr 'f' >> put_int_32 (rangeSize (bounds fa))
                   >> mapM_ (put_int_32 . fromFloat) (elems fa)
    IntArr   ia -> case put_some_int (elems ia) of
                    (c,op) -> putChr 'B' >> putChr c >> put_int_32 (rangeSize (bounds ia)-1)
                              >> mapM_ op (elems ia)
  where
    put_some_int :: [Int] -> (Char, Int -> Put)
    put_some_int is
        | all (between        0    0xff) is = ('C', put_int_8)
        | all (between   (-0x80)   0x7f) is = ('c', put_int_8)
        | all (between        0  0xffff) is = ('S', put_int_16)
        | all (between (-0x8000) 0x7fff) is = ('s', put_int_16)
        | all                      (> 0) is = ('I', put_int_32)
        | otherwise                         = ('i', put_int_32)

    between :: Int -> Int -> Int -> Bool
    between l r x = l <= x && x <= r

    fromFloat :: Float -> Int32
    fromFloat float = unsafePerformIO $ alloca $ \buf ->
                      poke (castPtr buf) float >> peek buf


isPaired, isProperlyPaired, isUnmapped, isMateUnmapped, isReversed,
    isMateReversed, isFirstMate, isSecondMate, isAuxillary, isFailsQC,
    isDuplicate, isTrimmed, isMerged :: BamRec -> Bool

isPaired         = flip testBit  0 . b_flag
isProperlyPaired = flip testBit  1 . b_flag
isUnmapped       = flip testBit  2 . b_flag
isMateUnmapped   = flip testBit  3 . b_flag
isReversed       = flip testBit  4 . b_flag
isMateReversed   = flip testBit  5 . b_flag
isFirstMate      = flip testBit  6 . b_flag
isSecondMate     = flip testBit  7 . b_flag
isAuxillary      = flip testBit  8 . b_flag
isFailsQC        = flip testBit  9 . b_flag
isDuplicate      = flip testBit 10 . b_flag
isTrimmed        = flip testBit 16 . b_flag
isMerged         = flip testBit 17 . b_flag


type_mask :: Int
type_mask = flagFirstMate .|. flagSecondMate .|. flagPaired


extAsInt :: Int -> String -> BamRec -> Int
extAsInt d nm br = case M.lookup nm (b_exts br) of Just (Int i) -> i ; _ -> d

extAsString :: String -> BamRec -> ByteString
extAsString nm br = case M.lookup nm (b_exts br) of
    Just (Char c) -> B.singleton c
    Just (Text s) -> s
    _             -> B.empty


setQualFlag :: Char -> BamRec -> BamRec
setQualFlag c br = br { b_exts = M.insert "ZQ" (Text s') $ b_exts br }
  where
    s  = extAsString "ZQ" br
    s' = if c `S.elem` s then s else c `S.cons` s


-- | Iteratee-style parser for SAM files, designed to be compatible with
-- the BAM parsers.  Parses plain uncompressed SAM, nothing else.  Since
-- it is supposed to work the same way as the BAM parser, it requires
-- the presense of the SQ header lines.  These are stripped from the
-- header text and turned into the symbol table.
decodeSam :: Monad m => (BamMeta -> Iteratee [BamRec] m a) -> Iteratee ByteString m (Iteratee [BamRec] m a)
decodeSam inner = joinI $ enumLinesBS $ do
    let pHeaderLine acc str = case P.parseOnly parseBamMetaLine str of Right f -> return $ f : acc
                                                                       Left e  -> fail $ e ++ ", " ++ show str
    meta <- liftM (foldr ($) mempty . reverse) (joinI $ I.breakE (not . S.isPrefixOf "@") $ I.foldM pHeaderLine [])
    decodeSamLoop (meta_refs meta) (inner meta)

decodeSamLoop :: Monad m => Refs -> Enumeratee [ByteString] [BamRec] m a
decodeSamLoop refs inner = I.convStream (liftI parse_record) inner
  where !refs' = M.fromList $ zip [ nm | BamSQ { sq_name = nm } <- F.toList refs ] [toEnum 0..]
        ref x = M.findWithDefault invalidRefseq x refs'

        parse_record (EOF x) = icont parse_record x
        parse_record (Chunk []) = liftI parse_record
        parse_record (Chunk (l:ls)) | "@" `S.isPrefixOf` l = parse_record (Chunk ls)
        parse_record (Chunk (l:ls)) = case P.parseOnly (parseSamRec ref) l of
            Right  r -> idone [fixup_bam_rec r] (Chunk ls)
            Left err -> icont parse_record (Just $ iterStrExc $ err ++ ", " ++ show l)

-- | Parser for SAM that doesn't look for a header.  Has the advantage
-- that it doesn't stall on a pipe that never delivers data.  Has the
-- disadvantage that it never reads the header and therefore needs a
-- list of allowed RNAMEs.
decodeSam' :: Monad m => Refs -> Enumeratee ByteString [BamRec] m a
decodeSam' refs inner = joinI $ enumLinesBS $ decodeSamLoop refs inner

parseSamRec :: (ByteString -> Refseq) -> P.Parser BamRec
parseSamRec ref = (\nm fl rn po mq cg rn' -> BamRec nm fl rn po mq cg (rn' rn))
                  <$> word <*> num <*> (ref <$> word) <*> (subtract 1 <$> num)
                  <*> (Q <$> num) <*> (Cigar <$> cigar) <*> rnext <*> (subtract 1 <$> num)
                  <*> snum <*> sequ <*> quals <*> exts <*> pure 0
  where
    sep      = P.endOfInput <|> () <$ P.char '\t'
    word     = P.takeTill ((==) '\t') <* sep
    num      = P.decimal <* sep
    snum     = P.signed P.decimal <* sep

    rnext    = id <$ P.char '=' <* sep <|> const . ref <$> word
    sequ     = {-# SCC "parseSamRec/sequ" #-}
               (V.empty <$ P.char '*' <|>
               V.fromList . map toNucleotides . S.unpack <$> P.takeWhile (P.inClass "acgtnACGTN")) <* sep

    quals    = {-# SCC "parseSamRec/quals" #-} B.empty <$ P.char '*' <* sep <|> B.map (subtract 33) <$> word

    cigar    = [] <$ P.char '*' <* sep <|>
               P.manyTill (flip (,) <$> P.decimal <*> cigop) sep

    cigop    = P.choice $ zipWith (\c r -> r <$ P.char c) "MIDNSHP" [Mat,Ins,Del,Nop,SMa,HMa,Pad]
    exts     = M.fromList <$> ext `P.sepBy` sep
    ext      = (\a b v -> ([a,b],v)) <$> P.anyChar <*> P.anyChar <*> (P.char ':' *> value)

    value    = P.char 'A' *> P.char ':' *> (Char <$>               anyWord8) <|>
               P.char 'i' *> P.char ':' *> (Int  <$>     P.signed P.decimal) <|>
               P.char 'Z' *> P.char ':' *> (Text <$> P.takeTill ((==) '\t')) <|>
               P.char 'H' *> P.char ':' *> (Bin  <$>               hexarray) <|>
               P.char 'f' *> P.char ':' *> (Float . realToFrac <$> P.double) <|>
               P.char 'B' *> P.char ':' *> (
                    P.satisfy (P.inClass "cCsSiI") *> (intArr   <$> many (P.char ',' *> P.signed P.decimal)) <|>
                    P.char 'f'                     *> (floatArr <$> many (P.char ',' *> P.double)))

    intArr   is = IntArr   $ listArray (0, length is -1) is
    floatArr fs = FloatArr $ listArray (0, length fs -1) $ map realToFrac fs
    hexarray    = B.pack . repack . S.unpack <$> P.takeWhile (P.inClass "0-9A-Fa-f")
    repack (a:b:cs) = fromIntegral (digitToInt a * 16 + digitToInt b) : repack cs ; repack _ = []

encodeSamEntry :: Refs -> BamRec -> String -> String
encodeSamEntry refs b = conjoin '\t' [
    unpck (b_qname b),
    shows (b_flag b .&. 0xffff),
    unpck (sq_name $ getRef refs $ b_rname b),
    shows (b_pos b + 1),
    shows (b_mapq b),
    shows (b_cigar b),
    unpck (sq_name $ getRef refs $ b_mrnm b),
    shows (b_mpos b + 1),
    shows (b_isize b + 1),
    shows (V.toList $ b_seq b),
    unpck (B.map (+33) $ b_qual b) ] .
    M.foldWithKey (\k v f -> (:) '\t' . (++) k . (:) ':' . extToSam v . f) id (b_exts b)
  where
    unpck = (++) . S.unpack
    conjoin c = foldr1 (\a f -> a . (:) c . f)

    extToSam (Int        i) = (:) 'i' . (:) ':' . shows i
    extToSam (Float      f) = (:) 'f' . (:) ':' . shows f
    extToSam (Text       t) = (:) 'Z' . (:) ':' . unpck t
    extToSam (Bin        x) = (:) 'H' . (:) ':' . tohex x
    extToSam (Char       c) = (:) 'A' . (:) ':' . (:) (w2c c)
    extToSam (IntArr   arr) = (:) 'B' . (:) ':' . (:) 'i' . sarr arr
    extToSam (FloatArr arr) = (:) 'B' . (:) ':' . (:) 'f' . sarr arr

    tohex = B.foldr (\c f -> w2d (c `shiftR` 4) . w2d (c .&. 0xf) . f) id
    w2d = (:) . S.index "0123456789ABCDEF" . fromIntegral
    sarr = conjoin ',' . map shows . elems
