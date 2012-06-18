{-# LANGUAGE OverloadedStrings, PatternGuards, BangPatterns, NoMonomorphismRestriction #-}

-- TODO:  
-- - Index writer
-- - Seeking, partially.  So far, we can only seek to the beginning of
--   the range of some RNAME.  Most of the functionality available from
--   the index is not yet supported.
-- - Writing of GZip'ed and plain BAM.  Probably more interesting as a
--   configurable wrapper.
-- - Automatic creation of some kind of index.  If possible, this should
--   be the standard index for sorted BAM.  Optionally a block index for
--   slicing of large files.  Maybe an index by name and an index for
--   group-sorted files.  All sensible indices should be generated
--   whenever a file is written.
-- - Same for statistics.  Something like "flagstats" could always be
--   written.  Actually, having @writeBamHandle@ return enhanced
--   flagstats as a result might be even better.
--
-- NICE_TO_HAVE:
-- - Decoding of BAM needlessly needs a MonadIO context, because of the
--   weird Zlib bindings.
-- - BZip compression isn't supported.  

-- TONOTDO:  
-- - SAM writer.  Writing SAM is a bad idea.
-- - Reader for gzipped/bzipped/bgzf'ed SAM.  Storing SAM is a bad idea,
--   so why would anyone ever want to compress, much less index it?

module Bio.File.Bam (
    module Bio.Base,

    Block,
    decompressBgzf,
    decompressBgzf',
    compressBgzf,

    BamrawEnumeratee,
    BamEnumeratee,

    isBam,
    isPlainBam,
    isGzipBam,
    isBgzfBam,
    isBamOrSam,

    decodeBam,
    decodeBamEntry,
    decodeBamSequence,

    decodeSam,
    decodeSam',

    decodeAnyBam,
    decodeAnyBamFile,
    decodeAnyBamOrSam,
    decodeAnyBamOrSamFile,

    encodeBam,
    encodeBamEntry,

    writeBamFile,
    writeBamHandle,

    BamRaw(..),

    BamRec(..),
    nullBamRec,
    Refs,
    noRefs,
    getRef,

    Refseq(..),
    invalidRefseq,
    isValidRefseq,
    invalidPos,
    isValidPos,
    unknownMapq,
    compareNames,

    MdOp(..),
    getMd,
    readMd,

    Cigar(..),
    CigOp(..),
    cigarToAlnLen,

    Extensions, Ext(..),
    extAsInt, extAsString, setQualFlag,

    flagPaired,         isPaired,
    flagProperlyPaired, isProperlyPaired,
    flagUnmapped,       isUnmapped,
    flagMateUnmapped,   isMateUnmapped,
    flagReversed,       isReversed,
    flagMateReversed,   isMateReversed,
    flagFirstMate,      isFirstMate,
    flagSecondMate,     isSecondMate,
    flagAuxillary,      isAuxillary,
    flagFailsQC,        isFailsQC,
    flagDuplicate,      isDuplicate,
    flagTrimmed,        isTrimmed,   
    flagMerged,         isMerged,       

    BamIndex,
    readBamIndex,
    readBamIndex',

    BamMeta(..),
    parseBamMeta,
    showBamMeta,

    BamHeader(..),
    BamSQ(..),
    BamSorting(..),
    BamOtherShit
) where

import Bio.Base
import Bio.File.Bgzf
import Bio.Iteratee

import Control.Monad
import Control.Applicative
import Data.Array.IArray
import Data.Array.IO
import Data.Array.Unboxed
import Data.Attoparsec              ( anyWord8 )
import Data.Attoparsec.Iteratee
import Data.Binary.Put
import Data.Bits                    ( testBit, shiftL, shiftR, (.&.), (.|.), complement )
import Data.Char                    ( chr, ord, isDigit, digitToInt )
import Data.Int                     ( Int64, Int32 )
import Data.Iteratee.ZLib
import Data.Monoid
import Data.Sequence                ( (<|), (|>), (><) )
import Data.Word                    ( Word32, Word8 )
import Foreign.Marshal.Alloc        ( alloca )
import Foreign.Ptr                  ( castPtr )
import Foreign.Storable             ( peek, poke )
import System.IO
import System.IO.Unsafe             ( unsafePerformIO )

import qualified Control.Monad.CatchIO          as CIO
import qualified Data.Attoparsec.Char8          as P
import qualified Data.Binary.Strict.Get         as G
import qualified Data.ByteString                as S
import qualified Data.ByteString.Char8          as B
import qualified Data.ByteString.Lazy.Char8     as L
import qualified Data.Foldable                  as F
import qualified Data.Iteratee                  as I
import qualified Data.Map                       as M
import qualified Data.Sequence                  as Z

-- ^ Parsers and Printers for BAM and SAM.  We employ an @Iteratee@
-- interface, and we strive to support everything possible in BAM.  So
-- far, the implementation of the nucleotides is somewhat lacking:  we
-- do not have support for ambiguity codes, and the "=" symbol is not
-- understood.


decompressBgzf' :: Monad m => Enumeratee S.ByteString Block m a
decompressBgzf' = decompress'

decompressBgzf :: Monad m => Enumeratee S.ByteString S.ByteString m a
decompressBgzf = decompress

compressBgzf :: Monad m => Enumeratee S.ByteString S.ByteString m a
compressBgzf = compress


-- | Cigar line in BAM coding
-- Bam encodes an operation and a length into a single integer, we keep
-- those integers in an array.
newtype Cigar = Cigar { unCigar :: [(CigOp, Int)] }

data CigOp = Mat | Ins | Del | Nop | SMa | HMa | Pad 
    deriving ( Eq, Ord, Enum, Show, Bounded, Ix )

instance Show Cigar where
    show (Cigar cs) = concat [ shows l (toChr op) | (op,l) <- cs ]
      where toChr = (:[]) . B.index "MIDNSHP" . fromEnum

-- | extracts the aligned length from a cigar line
-- This gives the length of an alignment as measured on the reference,
-- which is different from the length on the query or the length of the
-- alignment.
cigarToAlnLen :: Cigar -> Int
cigarToAlnLen (Cigar cig) = sum $ map l cig
  where l (op,n) = if op == Mat || op == Del || op == Nop then n else 0
    

-- | Reference sequence in Bam
-- Bam enumerates the reference sequences and then sorts by index.  We
-- need to track that index if we want to reproduce the sorting order.
newtype Refseq = Refseq { unRefseq :: Word32 } deriving (Show, Eq, Ord, Ix)

instance Enum Refseq where
    succ = Refseq . succ . unRefseq
    pred = Refseq . pred . unRefseq
    toEnum = Refseq . fromIntegral
    fromEnum = fromIntegral . unRefseq
    enumFrom = map Refseq . enumFrom . unRefseq
    enumFromThen (Refseq a) (Refseq b) = map Refseq $ enumFromThen a b
    enumFromTo (Refseq a) (Refseq b) = map Refseq $ enumFromTo a b
    enumFromThenTo (Refseq a) (Refseq b) (Refseq c) = map Refseq $ enumFromThenTo a b c


-- | Tests whether a reference sequence is valid.
-- Returns true unless the the argument equals @invalidRefseq@.
isValidRefseq :: Refseq -> Bool
isValidRefseq = (/=) invalidRefseq

-- | The invalid Refseq.
-- Bam uses this value to encode a missing reference sequence.
invalidRefseq :: Refseq
invalidRefseq = Refseq 0xffffffff

-- | The invalid Refseq.
-- Bam uses this value to encode a missing position.
invalidPos :: Int
invalidPos = -1

-- | Tests whether a position is valid.
-- Returns true unless the the argument equals @invalidPos@.
isValidPos :: Int -> Bool
isValidPos = (/=) invalidPos

unknownMapq :: Int
unknownMapq = 255

-- | internal representation of a BAM record
data BamRec = BamRec {
        b_qname :: S.ByteString,
        b_flag  :: Int,
        b_rname :: Refseq,
        b_pos   :: Int,
        b_mapq  :: Int,
        b_cigar :: Cigar,
        b_mrnm  :: Refseq,
        b_mpos  :: Int,
        b_isize :: Int,
        b_seq   :: [Nucleotide],
        b_qual  :: S.ByteString,       -- ^ quality, may be empty
        b_exts  :: Extensions,
        b_virtual_offset :: FileOffset -- ^ virtual offset for indexing purposes
    } deriving Show

nullBamRec :: BamRec
nullBamRec = BamRec {
        b_qname = S.empty,
        b_flag  = 0,
        b_rname = invalidRefseq,
        b_pos   = invalidPos,
        b_mapq  = unknownMapq,
        b_cigar = Cigar [],
        b_mrnm  = invalidRefseq,
        b_mpos  = invalidPos,
        b_isize = 0,
        b_seq   = [],
        b_qual  = S.empty,
        b_exts  = M.empty,
        b_virtual_offset = 0
    }

type BamrawEnumeratee m b = Enumeratee' BamMeta S.ByteString [BamRaw] m b
type BamEnumeratee m b = Enumeratee' BamMeta S.ByteString [BamRec] m b

-- | Tests if a data stream is a Bam file.
-- Recognizes plain Bam, gzipped Bam and bgzf'd Bam.  If a file is
-- recognized as Bam, a decoder (suitable Iteratee) for it is returned.
isBam, isPlainBam, isBgzfBam, isGzipBam :: MonadIO m => Iteratee S.ByteString m (Maybe (BamrawEnumeratee m a))
isBam = msum `liftM` sequence [ isPlainBam, isBgzfBam, isGzipBam ]

isPlainBam = (\n -> if n == 4 then Just ((>>= lift . run) . decompressPlain . decodeBam) else Nothing) 
             `liftM` i'lookAhead (I.heads "BAM\SOH")

isBgzfBam = (\n -> if n == 4 then Just ((>>= lift . run) . decompress' . decodeBam) else Nothing)
            `liftM` do b <- isBgzf
                       if b then i'lookAhead $ joinI $ decompress $ I.heads "BAM\SOH" else return 0

isGzipBam  = do b <- isGzip
                k <- if b then i'lookAhead $ joinI $ enumInflate GZip defaultDecompressParams isPlainBam
                          else return Nothing
                return $ (((>>= lift . run) . enumInflate GZip defaultDecompressParams) .) `fmap` k
                
isBamOrSam :: MonadIO m => Iteratee S.ByteString m (BamEnumeratee m a)
isBamOrSam = maybe decodeSam wrap `liftM` isBam
  where
    wrap enee it' = enee (\hdr -> I.mapStream decodeBamEntry (it' hdr)) >>= lift . run


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


-- | Checks if a file contains BAM in any of the common forms,
-- then decompresses it appropriately.  If the stream doesn't contain
-- BAM at all, it is instead decoded as SAM.  Since SAM is next to
-- impossible to recognize reliably, we don't even try.  Any old junk is
-- decoded as SAM and will fail later.
decodeAnyBamOrSam :: MonadIO m => BamEnumeratee m a
decodeAnyBamOrSam it = isBamOrSam >>= \k -> k it

decodeAnyBamOrSamFile :: MonadCatchIO m => FilePath -> (BamMeta -> Iteratee [BamRec] m a) -> m (Iteratee [BamRec] m a)
decodeAnyBamOrSamFile fn k = enumFileRandom defaultBufSize fn (decodeAnyBamOrSam k) >>= run


-- Seek to a given sequence in a Bam file, read those records.  This
-- requires an appropriate index (read separately), and the file must
-- have been opened in such a way as to allow seeking.  Enumerates over
-- the @BamRaw@ records of the correct sequence only, doesn't enumerate
-- at all if the sequence isn't found.

decodeBamSequence :: Monad m => BamIndex -> Refseq -> Enumeratee Block [BamRaw] m a
decodeBamSequence idx refseq iter = case idx ! refseq of
        _ | not (bounds idx `inRange` refseq) -> return iter
        0                                     -> return iter
        virtoff -> do virtualSeek $ fromIntegral virtoff
                      (decodeBamLoop ><> I.breakE wrong_ref) iter
  where
    wrong_ref br = let a = fromIntegral $ raw_data br `S.index` 0
                       b = fromIntegral $ raw_data br `S.index` 1
                       c = fromIntegral $ raw_data br `S.index` 2
                       d = fromIntegral $ raw_data br `S.index` 3
                       r = a `shiftL`  0 .|.  b `shiftL`  8 .|. 
                           c `shiftL` 16 .|.  d `shiftL` 24
                   in r /= unRefseq refseq
    

-- | A list of reference sequences.
type Refs = Z.Seq BamSQ

-- | The empty list of references.  Needed for BAM files that don't really store alignments.
noRefs :: Refs
noRefs = Z.empty

getRef :: Refs -> Refseq -> BamSQ
getRef refs (Refseq i) = Z.index refs (fromIntegral i)


-- | Decodes a raw block into a @BamRec@.
decodeBamEntry :: BamRaw -> BamRec
decodeBamEntry (BamRaw offs s) = case G.runGet go s of
    (Left  e,  _)             -> error e
    (Right r, s') | S.null s' -> fixup r
                  | otherwise -> error "incomplete BAM record"
  where
    go = do !rid       <- Refseq       <$> G.getWord32le
            !start     <- fromIntegral <$> G.getWord32le
            !namelen   <- fromIntegral <$> G.getWord8
            !mapq      <- fromIntegral <$> G.getWord8
            !_bin      <-                  G.getWord16le
            !cigar_len <- fromIntegral <$> G.getWord16le
            !flag      <- fromIntegral <$> G.getWord16le
            !read_len  <- fromIntegral <$> G.getWord32le
            !mate_rid  <- Refseq       <$> G.getWord32le
            !mate_pos  <- fromIntegral <$> G.getWord32le
            !ins_size  <- fromIntegral <$> G.getWord32le
            !read_name <- S.init       <$> G.getByteString namelen
            !cigar     <- Cigar . map decodeCigar <$> replicateM cigar_len G.getWord32le
            !qry_seq   <- G.getByteString $ (read_len+1) `div` 2
            !qual <- (\qs -> if S.all (0xff ==) qs then S.empty else qs) <$> G.getByteString read_len
            !exts <- getExtensions M.empty

            return $ BamRec read_name flag rid start mapq cigar
                            mate_rid mate_pos ins_size (take read_len $ expand qry_seq) qual exts offs
  
    bases = listArray (0,15) (map toNucleotide "NACNGNNNTNNNNNNN") :: Array Word8 Nucleotide
    expand t = if S.null t then [] else let x = S.head t in bases ! (x `shiftR` 4) : bases ! (x .&. 0xf) : expand (S.tail t)

    decodeCigar c | cc <= fromEnum (maxBound :: CigOp) = (toEnum cc, cl)
                  | otherwise = error "unknown Cigar operation"
      where cc = fromIntegral c .&. 0xf; cl = fromIntegral c `shiftR` 4

    -- fixups for changed conventions
    fixup b = (if b_flag b .&. flagLowQuality /= 0 then setQualFlag 'Q' else id) $
              (if b_flag b .&. flagLowComplexity /= 0 then setQualFlag 'C' else id) $
              b { b_flag = oflags .|. (eflags `shiftL` 16) }
      where
        flags' = b_flag b .&. complement (flagLowQuality .|. flagLowQuality)
        oflags | flags' .&. flagPaired == 0 = flags' .&. complement (flagFirstMate .|. flagSecondMate)
               | otherwise                  = flags'

        is_merged = flags' .&. (flagPaired .|. flagFirstMate .|. flagSecondMate) == (flagFirstMate .|. flagSecondMate)
        is_trimmed = flags' .&. (flagPaired .|. flagFirstMate .|. flagSecondMate) == flagSecondMate

        eflags = (if is_merged then flagMerged else 0) .|.
                 (if is_trimmed then flagTrimmed else 0) .|.
                 (case M.lookup "XF" (b_exts b) of Just (Int i) -> i ; _ -> 0)

        flagLowQuality = 0x800
        flagLowComplexity = 0x1000

-- | A collection of extension fields.  The key is actually only two @Char@s, but that proved impractical.
type Extensions = M.Map String Ext

data Ext = Int Int | Float Float | Text S.ByteString | Bin S.ByteString | Char Word8
         | IntArr (UArray Int Int) | FloatArr (UArray Int Float)
    deriving Show

getExtensions :: Extensions -> G.Get Extensions
getExtensions m = G.isEmpty >>= \e -> case e of
    True  -> return m
    False -> getExt >>= \(!k,!a) -> getExtensions $! M.insert k a m

getExt :: G.Get (String, Ext)
getExt = do key <- (\a b -> [w2c a, w2c b]) <$> G.getWord8 <*> G.getWord8
            typ <- G.getWord8
            res <- case w2c typ of
                    'Z' -> Text <$> getByteStringNul
                    'H' -> Bin  <$> getByteStringNul
                    'A' -> Char <$> G.getWord8
                    'f' -> Float . to_float <$> G.getWord32le
                    'B' -> get_arr
                    x | Just get <- M.lookup x get_some_int -> Int <$> get
                      | otherwise                           -> error $ "cannot handle optional field type " ++ [x]
            return (key,res)
  where
    to_float :: Word32 -> Float
    to_float word = unsafePerformIO $ alloca $ \buf ->
                    poke (castPtr buf) word >> peek buf

    get_arr = do tp <- chr . fromIntegral <$> G.getWord8
                 n <- fromIntegral <$> G.getWord32le
                 case tp of
                    'f' -> FloatArr . listArray (0,n) . map to_float <$> replicateM (n+1) G.getWord32le
                    _ | Just get <- M.lookup tp get_some_int -> IntArr . listArray (0,n) <$> replicateM (n+1) get
                      | otherwise                            -> error $ "cannot handle optional array field type " ++ [tp]

    get_some_int = M.fromList $ zip "cCsSiI" [
                        fromIntegral <$> G.getWord8,
                        fromIntegral <$> G.getWord8,
                        fromIntegral <$> G.getWord16le,
                        fromIntegral <$> G.getWord16le,
                        fromIntegral <$> G.getWord32le,
                        fromIntegral <$> G.getWord32le ]


        
getByteStringNul :: G.Get S.ByteString
getByteStringNul = S.init <$> (G.lookAhead (get_len 1) >>= G.getByteString)
  where 
    get_len l = G.getWord8 >>= \w -> if w == 0 then return l else get_len $! l+1


-- | Encode stuff into a BAM stream.
-- We send the encoded header and reference list to output through the
-- Bgzf compressor, then receive a list of records, which we concatenate
-- and send to output, too.
--
-- It would be nice if we were able to write an index on the side.  That
-- hasn't been designed in, yet.

encodeBam :: Monad m => BamMeta -> Enumeratee [S.ByteString] S.ByteString m a
encodeBam meta = eneeBam ><> compress
  where
    eneeBam = eneeCheckIfDone (\k -> eneeBam2 . k $ Chunk header)
    eneeBam2 = eneeCheckIfDone (\k -> eneeBam3 . k $ Chunk S.empty)
    eneeBam3 = eneeCheckIfDone (liftI . put)

    put k (EOF mx) = idone (liftI k) $ EOF mx
    put k (Chunk [    ]) = liftI $ put k
    put k (Chunk (r:rs)) = eneeCheckIfDone (\k' -> put k' (Chunk rs)) . k $ Chunk r'
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
                   put_int_32 . Z.length $ meta_refs meta
                   F.mapM_ putRef $ meta_refs meta

    putRef bs = do put_int_32 $ S.length (sq_name bs) + 1
                   putByteString $ sq_name bs
                   putWord8 0
                   put_int_32 $ sq_length bs


encodeBamEntry :: BamRec -> S.ByteString
encodeBamEntry = S.concat . L.toChunks . runPut . putEntry
  where
    putEntry  b = do putWord32le   $ unRefseq $ b_rname b
                     put_int_32    $ b_pos b
                     put_int_8     $ S.length (b_qname b) + 1
                     put_int_8     $ b_mapq b
                     put_int_16    $ distinctBin b
                     put_int_16    $ length $ unCigar $ b_cigar b
                     put_int_16    $ b_flag b
                     put_int_32    $ length $ b_seq b
                     putWord32le   $ unRefseq $ b_mrnm b
                     put_int_32    $ b_mpos b
                     put_int_32    $ b_isize b
                     putByteString $ b_qname b
                     putWord8 0
                     mapM_ (put_int_32 . encodeCigar) $ unCigar $ b_cigar b
                     putSeq $ b_seq b
                     putByteString $ if not (S.null (b_qual b)) then b_qual b
                                     else S.replicate (length $ b_seq b) 0xff
                     forM_ (M.toList $ more_exts b) $ \(k,v) -> 
                        case k of [c,d] -> putChr c >> putChr d >> putValue v
                                  _     -> error $ "invalid field key " ++ show k

    more_exts :: BamRec -> Extensions
    more_exts b = if xf /= 0 then x' else b_exts b
        where xf = b_flag b `shiftR` 16
              x' = M.insert "XF" (Int xf) $ b_exts b

    encodeCigar :: (CigOp,Int) -> Int
    encodeCigar (op,l) = fromEnum op .|. l `shiftL` 4

    putSeq :: [Nucleotide] -> Put
    putSeq (a:b:ns) = putWord8 (num a `shiftL` 4 .|. num b) >> putSeq ns
    putSeq [a]      = putWord8 (num a `shiftL` 4)
    putSeq [ ]      = return ()
  
    num :: Nucleotide -> Word8
    num A = 1 ; num C = 2 ; num G = 4 ; num T = 8 ; num N = 15 ; num _ = 0

-- | writes stuff to a BAM file                     
-- We generate BAM with dynamic blocks, then stream them out to the
-- file.
--
-- XXX This could write indexes on the side---a simple block index
-- for MapReduce style slicing, a standard BAM index or a name index
-- would be possible.
writeBamFile :: CIO.MonadCatchIO m => FilePath -> BamMeta -> Iteratee [BamRec] m ()
writeBamFile fp meta =
    CIO.bracket (liftIO $ openBinaryFile fp WriteMode)
                (liftIO . hClose)
                (flip writeBamHandle meta) 

writeBamHandle :: MonadIO m => Handle -> BamMeta -> Iteratee [BamRec] m ()
writeBamHandle hdl meta = 
    joinI $ I.mapStream encodeBamEntry $
    joinI $ encodeBam meta $
    mapChunksM_ (liftIO . S.hPut hdl)


put_int_32, put_int_16, put_int_8 :: Integral a => a -> Put
put_int_32 = putWord32le . fromIntegral
put_int_16 = putWord16le . fromIntegral
put_int_8  = putWord8 . fromIntegral

putChr :: Char -> Put
putChr = putWord8 . fromIntegral . ord

distinctBin :: BamRec -> Int
distinctBin b = mkbin 14 $ mkbin 17 $ mkbin 20 $ mkbin 23 $ mkbin 16 $ 0
  where beg = b_pos b
        end = beg + cigarToAlnLen (b_cigar b) - 1
        mkbin n x = if beg `shiftR` n /= end `shiftR` n then x
                    else ((1 `shiftL` (29-n))-1) `div` 7 + (beg `shiftR` n)

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


data MdOp = MdNum Int | MdRep Nucleotide | MdDel [Nucleotide] deriving Show

getMd :: BamRec -> Maybe [MdOp]
getMd r = case M.lookup "MD" $ b_exts r of
    Just (Text mdfield) -> readMd mdfield
    _                   -> Nothing

readMd :: B.ByteString -> Maybe [MdOp]
readMd s | B.null s           = return []
         | isDigit (B.head s) = do (n,t) <- B.readInt s
                                   (MdNum n :) <$> readMd t
         | B.head s == '^'    = let (a,b) = B.break isDigit (B.tail s)
                                in (MdDel (map toNucleotide $ B.unpack a) :) <$> readMd b
         | otherwise          = (MdRep (toNucleotide $ B.head s) :) <$> readMd (B.tail s)


flagPaired, flagProperlyPaired, flagUnmapped, flagMateUnmapped, flagReversed, flagMateReversed, flagFirstMate, flagSecondMate,
 flagAuxillary, flagFailsQC, flagDuplicate, flagTrimmed, flagMerged :: Int

flagPaired = 0x1
flagProperlyPaired = 0x2
flagUnmapped = 0x4
flagMateUnmapped = 0x8
flagReversed = 0x10
flagMateReversed = 0x20
flagFirstMate = 0x40
flagSecondMate = 0x80
flagAuxillary = 0x100
flagFailsQC = 0x200
flagDuplicate = 0x400

flagTrimmed = 0x10000
flagMerged  = 0x20000

isPaired, isProperlyPaired, isUnmapped, isMateUnmapped, isReversed,
    isMateReversed, isAuxillary, isFailsQC, isDuplicate, isTrimmed,
    isMerged :: BamRec -> Bool

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


-- | Compares two sequence names the way samtools does.
-- samtools sorts by "strnum_cmp":
-- . if both strings start with a digit, parse the initial
--   sequence of digits and compare numerically, if equal,
--   continue behind the numbers
-- . else compare the first characters (possibly NUL), if equal
--   continue behind them
-- . else both strings ended and the shorter one counts as
--   smaller (and that part is stupid)

compareNames :: Seqid -> Seqid -> Ordering
compareNames n m = case (B.uncons n, B.uncons m) of
        ( Nothing, Nothing ) -> EQ
        ( Just  _, Nothing ) -> GT
        ( Nothing, Just  _ ) -> LT
        ( Just (c,n'), Just (d,m') )
            | isDigit c && isDigit d -> 
                let Just (u,n'') = B.readInt n
                    Just (v,m'') = B.readInt m
                in case u `compare` v of 
                    LT -> LT
                    GT -> GT
                    EQ -> n'' `compareNames` m''
            | otherwise -> case c `compare` d of 
                    LT -> LT
                    GT -> GT
                    EQ -> n' `compareNames` m'
                                         

extAsInt :: Int -> String -> BamRec -> Int
extAsInt d nm br = case M.lookup nm (b_exts br) of Just (Int i) -> i ; _ -> d

extAsString :: String -> BamRec -> S.ByteString
extAsString nm br = case M.lookup nm (b_exts br) of
    Just (Char c) -> S.singleton c
    Just (Text s) -> s
    _             -> S.empty


setQualFlag :: Char -> BamRec -> BamRec
setQualFlag c br = br { b_exts = M.insert "ZQ" (Text s') $ b_exts br }
  where
    s  = extAsString "ZQ" br
    s' = if c `B.elem` s then s else c `B.cons` s

--
-- | Stop gap solution for a cheap index.  We only get the first offset
-- from the linear index, which allows us to navigate to a target
-- sequence.  Will do the rest when I need it.
type BamIndex = UArray Refseq Int64

readBamIndex :: FilePath -> IO BamIndex
readBamIndex = fileDriver readBamIndex'

readBamIndex' :: MonadIO m => Iteratee S.ByteString m BamIndex
readBamIndex' = do magic <- I.heads "BAI\1"
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

data BamMeta = BamMeta {
        meta_hdr :: BamHeader,
        meta_refs :: Refs,
        meta_other_shit :: [(Char, Char, BamOtherShit)],
        meta_comment :: [S.ByteString]
    } deriving Show

instance Monoid BamMeta where
    mempty = BamMeta mempty noRefs [] []
    a `mappend` b = BamMeta { meta_hdr = meta_hdr a `mappend` meta_hdr b
                            , meta_refs = meta_refs a >< meta_refs b
                            , meta_other_shit = meta_other_shit a ++ meta_other_shit b
                            , meta_comment = meta_comment a ++ meta_comment b }

data BamHeader = BamHeader {
        hdr_version :: (Int, Int),
        hdr_sorting :: BamSorting,
        hdr_other_shit :: BamOtherShit
    } deriving Show

instance Monoid BamHeader where
    mempty = BamHeader (0,0) Unsorted []
    a `mappend` b = BamHeader { hdr_version = hdr_version a `min` hdr_version b
                              , hdr_sorting = let u = hdr_sorting a ; v = hdr_sorting b in if u == v then u else Unsorted
                              , hdr_other_shit = hdr_other_shit a ++ hdr_other_shit b }

data BamSQ = BamSQ {
        sq_name :: Seqid,
        sq_length :: Int,
        sq_other_shit :: BamOtherShit
    } deriving Show

bad_seq :: BamSQ
bad_seq = BamSQ (error "no SN field") (error "no LN field") []

data BamSorting = Unsorted | Grouped | Queryname | Coordinate | GroupSorted 
    deriving (Show, Eq)

type BamOtherShit = [(Char, Char, S.ByteString)]

parseBamMeta :: P.Parser BamMeta
parseBamMeta = foldr ($) mempty <$> P.sepBy parseBamMetaLine (P.char '\n')

parseBamMetaLine :: P.Parser (BamMeta -> BamMeta)
parseBamMetaLine = P.char '@' >> P.choice [hdLine, sqLine, coLine, otherLine]
  where
    hdLine = P.string "HD\t" >> 
             (\fns meta -> meta { meta_hdr = foldr ($) (meta_hdr meta) fns })
               <$> P.sepBy1 (P.choice [hdvn, hdso, hdother]) (P.char '\t')
    
    sqLine = P.string "SQ\t" >> 
             (\fns meta -> meta { meta_refs = foldr ($) bad_seq fns <| meta_refs meta })
               <$> P.sepBy1 (P.choice [sqnm, sqln, sqother]) (P.char '\t')
    
    hdvn = P.string "VN:" >>
           (\a b hdr -> hdr { hdr_version = (a,b) })
             <$> P.decimal <*> ((P.char '.' <|> P.char ':') >> P.decimal)

    hdso = P.string "SO:" >>
           (\s hdr -> hdr { hdr_sorting = s })
             <$> P.choice [ Grouped  <$ P.string "grouped"
                          , Queryname <$ P.string "queryname"
                          , Coordinate <$ P.string "coordinate"
                          , GroupSorted <$ P.string "groupsort"
                          , Unsorted <$ P.skipWhile (\c -> c/='\t' && c/='\n') ]

    sqnm = P.string "SN:" >> (\s sq -> sq { sq_name = s }) <$> pall
    sqln = P.string "LN:" >> (\i sq -> sq { sq_length = i }) <$> P.decimal

    hdother = (\t hdr -> hdr { hdr_other_shit = t : hdr_other_shit hdr }) <$> tagother
    sqother = (\t sq  -> sq  { sq_other_shit = t : sq_other_shit sq }) <$> tagother
    
    coLine = P.string "CO\t" >>
             (\s meta -> meta { meta_comment = s : meta_comment meta })
               <$> P.takeWhile (/= 'n')

    otherLine = (\a b ts meta -> meta { meta_other_shit = (a,b,ts) : meta_other_shit meta })
                  <$> P.anyChar <*> P.anyChar <*> (P.char '\t' >> P.sepBy1 tagother (P.char '\t'))

    tagother :: P.Parser (Char,Char,S.ByteString)
    tagother = (,,) <$> P.anyChar <*> P.anyChar <*> (P.char ':' >> pall)
    
    pall :: P.Parser S.ByteString
    pall = P.takeWhile (\c -> c/='\t' && c/='\n')

showBamMeta :: BamMeta -> L.ByteString -> L.ByteString
showBamMeta (BamMeta h ss os cs) = 
    show_bam_meta_hdr h .
    F.foldr ((.) . show_bam_meta_seq) id ss .
    foldr ((.) . show_bam_meta_other) id os .
    foldr ((.) . show_bam_meta_comment) id cs
  where
    show_bam_meta_hdr (BamHeader (major,minor) so os') = 
        L.append "@HD\tVN:" . L.append (L.pack (show major ++ '.' : show minor)) .
        L.append (case so of Unsorted -> L.empty
                             Grouped  -> "\tSO:grouped"
                             Queryname  -> "\tSO:queryname"
                             Coordinate  -> "\tSO:coordinate"
                             GroupSorted  -> "\tSO:groupsort") .
        show_bam_others os'

    show_bam_meta_seq (BamSQ  _  _ []) = id
    show_bam_meta_seq (BamSQ nm ln ts) =
        L.append "@SQ\tSN:" . L.append (L.fromChunks [nm]) . L.append "\tLN:" .
        L.append (L.pack (show ln)) . show_bam_others ts

    show_bam_meta_comment cm = L.append "@CO\t" . L.append (L.fromChunks [cm]) . L.cons '\n'

    show_bam_meta_other (a,b,ts) = 
        L.cons '@' . L.cons a . L.cons b . show_bam_others ts

    show_bam_others ts =         
        foldr ((.) . show_bam_other) id ts . L.cons '\n'

    show_bam_other (a,b,v) = 
        L.cons '\t' . L.cons a . L.cons b . L.cons ':' . L.append (L.fromChunks [v])


-- | Bam record in its native encoding along with virtual address.
data BamRaw = BamRaw { virt_offset :: {-# UNPACK #-} !FileOffset
                     , raw_data :: {-# UNPACK #-} !S.ByteString }


-- | Decode a BAM stream into raw entries.  Note that the entries can be
-- unpacked using @decodeBamEntry@.  Also note that this is an
-- Enumeratee in spirit, only the @BamMeta@ and @Refs@ need to get
-- passed separately.
decodeBam :: Monad m => (BamMeta -> Iteratee [BamRaw] m a) -> Iteratee Block m (Iteratee [BamRaw] m a)
decodeBam inner = do meta <- liftBlock get_bam_header
                     refs <- liftBlock get_ref_array
                     decodeBamLoop $ inner $! merge meta refs
  where
    get_bam_header  = do magic <- I.heads "BAM\SOH"
                         when (magic /= 4) $ fail "BAM signature not found"
                         hdr_len <- endianRead4 LSB
                         joinI $ I.take (fromIntegral hdr_len) $ parserToIteratee parseBamMeta

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
    loop k = I.isFinished >>= loop' k
    loop' k True = return $ liftI k
    loop' k False = do off <- getOffset
                       raw <- liftBlock $ do
                                bsize <- endianRead4 LSB
                                when (bsize < 32) $ fail "short BAM record"
                                i'getString (fromIntegral bsize)
                       eneeCheckIfDone loop . k $ Chunk [BamRaw off raw]


-- | Iteratee-style parser for SAM files, designed to be compatible with
-- the BAM parsers.  Parses plain uncompressed SAM, nothing else.  Since
-- it is supposed to work the same way as the BAM parser, it requires
-- the presense of the SQ header lines.  These are stripped from the
-- header text and turned into the symbol table.
decodeSam :: Monad m => (BamMeta -> Iteratee [BamRec] m a) -> Iteratee S.ByteString m (Iteratee [BamRec] m a)
decodeSam inner = joinI $ enumLinesBS $ do
    let pHeaderLine acc str = case P.parseOnly parseBamMetaLine str of Right f -> return $ f : acc
                                                                       Left e  -> fail $ e ++ ", " ++ show str
    meta <- liftM (foldr ($) mempty . reverse) (joinI $ I.breakE (not . S.isPrefixOf "@") $ I.foldM pHeaderLine [])
    decodeSamLoop (meta_refs meta) (inner meta)


decodeSamLoop :: Monad m => Refs -> Enumeratee [S.ByteString] [BamRec] m a
decodeSamLoop refs inner = I.convStream (liftI parse_record) inner
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
decodeSam' :: Monad m => Refs -> Enumeratee S.ByteString [BamRec] m a
decodeSam' refs inner = joinI $ enumLinesBS $ decodeSamLoop refs inner

parseSamRec :: (B.ByteString -> Refseq) -> P.Parser BamRec
parseSamRec ref = (\nm fl rn po mq cg rn' -> BamRec nm fl rn po mq cg (rn' rn))
                  <$> word <*> num <*> (ref <$> word) <*> (subtract 1 <$> num)
                  <*> num <*> (Cigar <$> cigar) <*> rnext <*> (subtract 1 <$> num)
                  <*> snum <*> sequ <*> quals <*> exts <*> pure 0
  where
    sep      = P.endOfInput <|> () <$ P.char '\t'
    word     = P.takeTill ((==) '\t') <* sep
    num      = P.decimal <* sep
    snum     = P.signed P.decimal <* sep

    rnext    = id <$ P.char '=' <* sep <|> const . ref <$> word
    sequ     = {-# SCC "parseSamRec/sequ" #-}
               ([] <$ P.char '*' <|>
               map toNucleotide . B.unpack <$> P.takeWhile (P.inClass "acgtnACGTN")) <* sep
    
    quals    = {-# SCC "parseSamRec/quals" #-} S.empty <$ P.char '*' <* sep <|> S.map (subtract 33) <$> word

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
    hexarray    = S.pack . repack . B.unpack <$> P.takeWhile (P.inClass "0-9A-Fa-f")
    repack (a:b:cs) = fromIntegral (digitToInt a * 16 + digitToInt b) : repack cs ; repack _ = []

-- ------------------------------------------------------------------- Tests

some_file :: FilePath
some_file = "/mnt/ngs_data/101203_SOLEXA-GA04_00007_PEDi_MM_QF_SR/Ibis/BWA/s_5_L3280_sequence_mq_hg19_nohap.bam"

bam_test' :: FilePath -> IO ()
bam_test' = fileDriver $
            joinI $ decompress'  $
            joinI $ decodeBam    $
            dump_bam
            
bam_test :: FilePath -> IO ()
bam_test = fileDriverRandom $
           joinI $ decompress'  $
           joinI $ do virtualSeek 0
                      decodeBam dump_bam 

dump_bam :: BamMeta -> Iteratee [BamRaw] IO ()
dump_bam meta = lift (print meta) >> print_names

seek_test :: [Char] -> Word32 -> IO ()
seek_test fp i = do
    idx <- readBamIndex $ fp ++ ".bai"
    flip fileDriverRandom fp $
           joinI $ decompress'  $
           joinI $ decodeBamSequence idx (Refseq i) print_names_and_refs

sam_test :: IO ()
sam_test = fileDriver (joinI $ decodeSam (const print_names')) "foo.sam"

print_names :: Iteratee [BamRaw] IO ()
print_names = I.mapM_ $ B.putStrLn . b_qname . decodeBamEntry 

print_names_and_refs :: Iteratee [BamRaw] IO ()
print_names_and_refs = I.mapM_ $ pr . decodeBamEntry
  where pr b = putStrLn $ shows (b_qname b) " " ++ show (b_rname b)

print_names' :: Iteratee [BamRec] IO ()
print_names' = I.mapM_ $ B.putStrLn . b_qname


bam2bam_test :: IO ()
bam2bam_test = withFile "foo.bam" WriteMode $       \hdl ->
               flip fileDriver some_file $
               joinI $ decompress' $ 
               joinI $ decodeBam   $                \meta ->
               joinI $ I.mapStream raw_data $
               joinI $ encodeBam meta $
               mapChunksM_ (S.hPut hdl)

sam2bam_test :: IO ()               
sam2bam_test = withFile "bar.bam" WriteMode $                   \hdl ->
               flip fileDriver "foo.sam" $
               joinI $ decodeSam $                              \meta ->
               joinI $ I.mapStream encodeBamEntry $
               lift (print meta) >>=                            \_ -> 
               joinI $ encodeBam meta $
               mapChunksM_ (S.hPut hdl)


