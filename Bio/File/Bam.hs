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
-- - Reader for gzipped/bzipped/bgzf'ed SAM.  Storing SAM is a bad idea,
--   so why would anyone ever want to compress, much less index it?

module Bio.File.Bam (
    module Bio.Base,
    module Bio.File.Bam.Header,
    module Bio.File.Bam.Raw,

    Block(..),
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

    MdOp(..),
    getMd,
    readMd,
    showMd,

    Cigar(..),
    CigOp(..),
    cigarToAlnLen,
    cigar_op,
    cigar_len,

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

    Word32, Word8
) where

import Bio.Base
import Bio.File.Bam.Header
import Bio.File.Bam.Raw
import Bio.Iteratee

import Control.Monad
import Control.Applicative
import Data.Array.IArray
import Data.Array.Unboxed
import Data.Attoparsec              ( anyWord8 )
import Data.Binary.Put
import Data.Bits                    ( Bits, testBit, shiftL, shiftR, (.&.), (.|.), complement )
import Data.Char                    ( ord, isDigit, digitToInt )
import Data.Int                     ( Int32 )
import Data.Monoid                  ( mempty )
import Data.Vector.Generic          ( (!?) )
import Data.Word                    ( Word32, Word8 )
import Foreign.Marshal.Alloc        ( alloca )
import Foreign.Ptr                  ( castPtr )
import Foreign.Storable             ( peek, poke )
import System.IO
import System.IO.Unsafe             ( unsafePerformIO )

import qualified Data.Attoparsec.Char8          as P
import qualified Data.Binary.Strict.Get         as G
import qualified Data.ByteString                as B
import qualified Data.ByteString.Char8          as S
import qualified Data.ByteString.Lazy.Char8     as L
import qualified Data.Foldable                  as F
import qualified Data.Iteratee                  as I
import qualified Data.Map                       as M
import qualified Data.Vector.Generic            as V

-- ^ Parsers and Printers for BAM and SAM.  We employ an @Iteratee@
-- interface, and we strive to support everything possible in BAM.  So
-- far, the implementation of the nucleotides is somewhat lacking:  we
-- do not have support for ambiguity codes, and the "=" symbol is not
-- understood.

type ByteString = B.ByteString

-- | Cigar line in BAM coding
-- Bam encodes an operation and a length into a single integer, we keep
-- those integers in an array.
newtype Cigar = Cigar { unCigar :: [(CigOp, Int)] }

data CigOp = Mat | Ins | Del | Nop | SMa | HMa | Pad
    deriving ( Eq, Ord, Enum, Show, Bounded, Ix )

instance Show Cigar where
    show (Cigar cs) = concat [ shows l (toChr op) | (op,l) <- cs ]
      where toChr = (:[]) . S.index "MIDNSHP=X" . fromEnum

-- | extracts the aligned length from a cigar line
-- This gives the length of an alignment as measured on the reference,
-- which is different from the length on the query or the length of the
-- alignment.
cigarToAlnLen :: Cigar -> Int
cigarToAlnLen (Cigar cig) = sum $ map l cig
  where l (op,n) = if op == Mat || op == Del || op == Nop then n else 0
    
cigar_op :: Word32 -> CigOp
cigar_op x = toEnum $ fromIntegral $ x .&. 0xf

cigar_len :: Word32 -> Int
cigar_len x = fromIntegral x `shiftR` 4

-- | internal representation of a BAM record
data BamRec = BamRec {
        b_qname :: !Seqid,
        b_flag  :: !Int,
        b_rname :: !Refseq,
        b_pos   :: !Int,
        b_mapq  :: !Int,
        b_cigar :: !Cigar,
        b_mrnm  :: !Refseq,
        b_mpos  :: !Int,
        b_isize :: !Int,
        b_seq   :: !Sequence,
        b_qual  :: !ByteString,         -- ^ quality, may be empty
        b_exts  :: Extensions,
        b_virtual_offset :: !FileOffset -- ^ virtual offset for indexing purposes
    } deriving Show

nullBamRec :: BamRec
nullBamRec = BamRec {
        b_qname = S.empty,
        b_flag  = 0,
        b_rname = invalidRefseq,
        b_pos   = invalidPos,
        b_mapq  = 0,
        b_cigar = Cigar [],
        b_mrnm  = invalidRefseq,
        b_mpos  = invalidPos,
        b_isize = 0,
        b_seq   = V.empty,
        b_qual  = S.empty,
        b_exts  = M.empty,
        b_virtual_offset = 0
    }

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

decodeAnyBamOrSamFile :: MonadCatchIO m => FilePath -> (BamMeta -> Iteratee [BamRec] m a) -> m (Iteratee [BamRec] m a)
decodeAnyBamOrSamFile fn k = enumFileRandom defaultBufSize fn (decodeAnyBamOrSam k) >>= run


    


-- | Decodes a raw block into a @BamRec@.
decodeBamEntry :: BamRaw -> BamRec
decodeBamEntry br = either error fixup . fst . G.runGet go $ raw_data br
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
            !qual <- (\qs -> if B.all (0xff ==) qs then B.empty else qs) <$> G.getByteString read_len
            !exts <- getExtensions M.empty

            return $ BamRec read_name flag rid start mapq cigar
                            mate_rid mate_pos ins_size (V.fromListN read_len $ expand qry_seq) qual exts (virt_offset br)
  
    bases = listArray (0,15) (map toNucleotide "NACNGNNNTNNNNNNN") :: Array Word8 Nucleotide
    expand t = if S.null t then [] else let x = B.head t in bases ! (x `shiftR` 4) : bases ! (x .&. 0xf) : expand (B.tail t)

    decodeCigar c | cc <= fromEnum (maxBound :: CigOp) = (toEnum cc, cl)
                  | otherwise = error "unknown Cigar operation"
      where cc = fromIntegral c .&. 0xf; cl = fromIntegral c `shiftR` 4

    -- fixups for changed conventions
    fixup b = (if b_flag b .&. flagLowQuality /= 0 then setQualFlag 'Q' else id) $          -- low qual, new convention
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

getExtensions :: Extensions -> G.Get Extensions
getExtensions m = getExt `G.plus` return m
  where
    getExt :: G.Get Extensions
    getExt = do
            key <- (\a b -> [w2c a, w2c b]) <$> G.getWord8 <*> G.getWord8
            typ <- G.getWord8
            let cont v = getExtensions $! M.insert key v m
            case w2c typ of
                    'Z' -> cont . Text =<< getByteStringNul
                    'H' -> cont . Bin  =<< getByteStringNul
                    'A' -> cont . Char =<< G.getWord8
                    'f' -> cont . Float . to_float =<< G.getWord32le
                    'B' -> do tp <- G.getWord8
                              n <- fromIntegral <$> G.getWord32le
                              case w2c tp of
                                 'f' -> cont . FloatArr . listArray (0,n) . map to_float =<< replicateM (n+1) G.getWord32le
                                 x | Just get <- M.lookup x get_some_int -> cont . IntArr . listArray (0,n) =<< replicateM (n+1) get
                                   | otherwise                           -> fail $ "array type code " ++ show x ++ " not recognized"
                    x | Just get <- M.lookup x get_some_int -> cont . Int =<< get
                      | otherwise                           -> fail $ "type code " ++ show x ++ " not recognized"

    to_float :: Word32 -> Float
    to_float word = unsafePerformIO $ alloca $ \buf ->
                    poke (castPtr buf) word >> peek buf

    get_some_int :: M.Map Char (G.Get Int)
    get_some_int = M.fromList $ zip "cCsSiI" [
                        fromIntegral <$> G.getWord8,
                        fromIntegral <$> G.getWord8,
                        fromIntegral <$> G.getWord16le,
                        fromIntegral <$> G.getWord16le,
                        fromIntegral <$> G.getWord32le,
                        fromIntegral <$> G.getWord32le ]


        
getByteStringNul :: G.Get ByteString
getByteStringNul = S.init <$> (G.lookAhead (get_len 1) >>= G.getByteString)
  where 
    get_len l = G.getWord8 >>= \w -> if w == 0 then return l else get_len $! l+1




encodeBamEntry :: BamRec -> BamRaw
encodeBamEntry = bamRaw 0 . S.concat . L.toChunks . runPut . putEntry
  where
    putEntry  b = do putWord32le   $ unRefseq $ b_rname b
                     put_int_32    $ b_pos b
                     put_int_8     $ S.length (b_qname b) + 1
                     put_int_8     $ b_mapq b
                     put_int_16    $ distinctBin b
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

    putSeq :: Sequence -> Put
    putSeq v = case v !? 0 of 
                 Nothing -> return ()
                 Just a  -> case v !? 1 of
                    Nothing -> putWord8 (num a `shiftL` 4)
                    Just b  -> putWord8 (unN a `shiftL` 4 .|. num b)
                               >> putSeq (V.drop 2 v)
                         
    num :: Nucleotide -> Word8
    num (N x) = x .&. 15

-- | writes BAM encoded stuff to a @Handle@
-- We generate BAM with dynamic blocks, then stream them out to the file.
--
-- XXX This could write indexes on the side---a simple block index
-- for MapReduce style slicing, a standard BAM index or a name index
-- would be possible.
writeBamHandle :: MonadIO m => Handle -> BamMeta -> Iteratee [BamRec] m ()
writeBamHandle hdl meta = I.mapStream encodeBamEntry =$ writeRawBamHandle hdl meta

-- | writes BAM encoded stuff to a file
-- XXX This should(!) write indexes on the side---a simple block index
-- for MapReduce style slicing, a standard BAM index or a name index
-- would be possible.  When writing to a file, this makes even more
-- sense than when writing to a @Handle@.
writeBamFile :: MonadCatchIO m => FilePath -> BamMeta -> Iteratee [BamRec] m ()
writeBamFile fp meta = I.mapStream encodeBamEntry =$ writeRawBamFile fp meta

-- | write BAM encoded stuff to stdout
-- This send uncompressed BAM to stdout.  Useful for piping to other
-- tools.
pipeBamOutput :: MonadIO m => BamMeta -> Iteratee [BamRec] m ()
pipeBamOutput meta = I.mapStream encodeBamEntry =$ pipeRawBamOutput meta

-- | write in SAM format to stdout
-- This is useful for piping to other tools (say, AWK scripts) or for
-- debugging.  No convenience function to send SAM to a file exists,
-- because that's a stupid idea.
pipeSamOutput :: MonadIO m => BamMeta -> Iteratee [BamRec] m ()
pipeSamOutput meta = do liftIO . L.putStr $ showBamMeta meta L.empty
                        mapStreamM_ $ \b -> liftIO . putStr $ encodeSamEntry (meta_refs meta) b "\n"

pipeRawSamOutput :: MonadIO m => BamMeta -> Iteratee [BamRaw] m ()
pipeRawSamOutput hdr = joinI $ mapStream decodeBamEntry $ pipeSamOutput hdr

put_int_32, put_int_16, put_int_8 :: Integral a => a -> Put
put_int_32 = putWord32le . fromIntegral
put_int_16 = putWord16le . fromIntegral
put_int_8  = putWord8 . fromIntegral

putChr :: Char -> Put
putChr = putWord8 . fromIntegral . ord

distinctBin :: BamRec -> Int
distinctBin b = mkbin 14 $ mkbin 17 $ mkbin 20 $ mkbin 23 $ mkbin 26 $ 0
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
    Just (Char mdfield) -> readMd $ B.singleton mdfield
    _                   -> Nothing

readMd :: ByteString -> Maybe [MdOp]
readMd s | S.null s           = return []
         | isDigit (S.head s) = do (n,t) <- S.readInt s
                                   (MdNum n :) <$> readMd t
         | S.head s == '^'    = let (a,b) = S.break isDigit (S.tail s)
                                in (MdDel (map toNucleotide $ S.unpack a) :) <$> readMd b
         | otherwise          = (MdRep (toNucleotide $ S.head s) :) <$> readMd (S.tail s)

showMd :: [MdOp] -> ByteString
showMd = S.pack . flip s1 []
  where
    s1 (MdNum  i : MdNum  j : ms) = s1 (MdNum (i+j) : ms)
    s1 (MdNum  0            : ms) = s1 ms
    s1 (MdNum  i            : ms) = shows i . s1 ms

    s1 (MdRep  r            : ms) = shows r . s1 ms

    s1 (MdDel d1 : MdDel d2 : ms) = s1 (MdDel (d1++d2) : ms)
    s1 (MdDel []            : ms) = s1 ms
    s1 (MdDel ns : MdRep  r : ms) = (:) '^' . shows ns . (:) '0' . shows r . s1 ms
    s1 (MdDel ns            : ms) = (:) '^' . shows ns . s1 ms
    s1 [                        ] = id


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
            Right  r -> idone [r] (Chunk ls)
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
                  <*> num <*> (Cigar <$> cigar) <*> rnext <*> (subtract 1 <$> num)
                  <*> snum <*> sequ <*> quals <*> exts <*> pure 0
  where
    sep      = P.endOfInput <|> () <$ P.char '\t'
    word     = P.takeTill ((==) '\t') <* sep
    num      = P.decimal <* sep
    snum     = P.signed P.decimal <* sep

    rnext    = id <$ P.char '=' <* sep <|> const . ref <$> word
    sequ     = {-# SCC "parseSamRec/sequ" #-}
               (V.empty <$ P.char '*' <|>
               V.fromList . map toNucleotide . S.unpack <$> P.takeWhile (P.inClass "acgtnACGTN")) <* sep
    
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

-- ------------------------------------------------------------------- Tests

some_file :: FilePath
some_file = "/mnt/ngs_data/101203_SOLEXA-GA04_00007_PEDi_MM_QF_SR/Ibis/BWA/s_5_L3280_sequence_mq_hg19_nohap.bam"

bam_test' :: FilePath -> IO ()
bam_test' = fileDriver $
            joinI $ decompressBgzf $
            joinI $ decodeBam      $
            dump_bam
            
bam_test :: FilePath -> IO ()
bam_test = fileDriverRandom $
           joinI $ decompressBgzf $
           joinI $ seek 0 >> decodeBam dump_bam

dump_bam :: BamMeta -> Iteratee [BamRaw] IO ()
dump_bam meta = lift (print meta) >> print_names

seek_test :: [Char] -> Word32 -> IO ()
seek_test fp i = do
    idx <- readBamIndex $ fp ++ ".bai"
    flip fileDriverRandom fp $
           joinI $ decompressBgzf $
           joinI $ decodeBamSequence idx (Refseq i) print_names_and_refs

sam_test :: IO ()
sam_test = fileDriver (joinI $ decodeSam (const print_names')) "foo.sam"

print_names :: Iteratee [BamRaw] IO ()
print_names = I.mapM_ $ S.putStrLn . br_qname

print_names_and_refs :: Iteratee [BamRaw] IO ()
print_names_and_refs = I.mapM_ pr
  where pr b = putStrLn $ shows (br_qname b) " " ++ show (br_rname b)

print_names' :: Iteratee [BamRec] IO ()
print_names' = I.mapM_ $ S.putStrLn . b_qname


bam2bam_test :: IO ()
bam2bam_test = withFile "foo.bam" WriteMode $       \hdl ->
               flip fileDriver some_file    $
               joinI $ decompressBgzf       $
               joinI $ decodeBam            $       \meta ->
               joinI $ encodeBam meta       $
               mapChunksM_ (S.hPut hdl)

sam2bam_test :: IO ()               
sam2bam_test = withFile "bar.bam" WriteMode       $             \hdl ->
               flip fileDriver "foo.sam"          $
               joinI $ decodeSam                  $             \meta ->
               joinI $ I.mapStream encodeBamEntry $
               lift (print meta)                >>=             \_ -> 
               joinI $ encodeBam meta             $
               mapChunksM_ (S.hPut hdl)


