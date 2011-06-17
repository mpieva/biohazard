{-# LANGUAGE OverloadedStrings, PatternGuards, BangPatterns #-}
module Bio.File.Bam (
    module Bio.Base,

    -- isBam,

    decodeBam,
    decodeBamEntry,
    encodeBam,
    encodeBamEntry,

    -- writeBamFile,

    BamRaw(..),

    BamRec(..),
    Refs,
    noRefs,
    (!),

    Refseq(..),
    invalidRefseq,
    isValidRefseq,
    invalidPos,
    isValidPos,
    compareNames,

    MdOp(..),
    getMd,
    readMd,

    Cigar(..),
    CigOp(..),
    cigarToAlnLen,

    Extensions,
    Ext(..),

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
    flagLowQuality,     isLowQuality,
    flagLowComplexity,  isLowComplexity,

    filterFlags,
    isMerged,
    isAdapterTrimmed,

    -- BamIndex,
    -- readBamIndex,
    -- readBamSequence,

    BamMeta(..),
    nullMeta,
    parseBamMeta,
    showBamMeta,

    BamHeader(..),
    BamSorting(..),
    BamOtherShit
) where

import Bio.Base
import Bio.File.Bgzf

import Control.Monad                ( replicateM, replicateM_, when, forM_ )
import Control.Monad.IO.Class
import Control.Monad.Trans.Class
import Control.Applicative          ( (<$>), (<$), (<*>), (<*) )
import Data.Array.IO
import Data.Array.Unboxed
import Data.Attoparsec.Iteratee
import Data.Binary.Get
import Data.Binary.Put
import Data.Bits                    ( testBit, shiftL, shiftR, (.&.), (.|.) )
import Data.Char                    ( chr, ord, isDigit )
import Data.Int                     ( Int64, Int32 )
import Data.Iteratee.Base
import Data.Iteratee.Binary
import Data.Word                    ( Word64, Word32, Word8 )
import Foreign.Marshal.Alloc        ( alloca )
import Foreign.Ptr                  ( castPtr )
import Foreign.Storable             ( peek, poke )
import System.IO.Unsafe             ( unsafePerformIO )

import qualified Data.Attoparsec.Char8          as P
import qualified Data.Binary.Strict.Get         as G
import qualified Data.ByteString                as S
import qualified Data.ByteString.Char8          as B
import qualified Data.ByteString.Lazy.Char8     as L
import qualified Data.Iteratee                  as I
import qualified Data.Map                       as M


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
        b_virtual_offset :: Int64      -- ^ virtual offset for indexing purposes
    } deriving Show

{-
-- | Tests if a data stream is a Bam file.
-- Recognizes plain Bam, gzipped Bam and bgzf'd Bam.
isBam :: L.ByteString -> Bool
isBam s = L.pack "BAM\SOH" `L.isPrefixOf` s ||
           ( isGzip s && isBam (decompressBgzf s))

-- | Tests if a stream is a GZip stream.
-- This only tests the magic number.  Since Bgzf (and therefore Bam) is
-- GZip with added conventions, these files also return true.
isGzip :: L.ByteString -> Bool
isGzip s = not (L.null (L.drop 26 s)) && L.pack "\31\139" `L.isPrefixOf` s
-}


{-
-- Seek to a given sequence in a Bam file, read records.  Opens a new
-- handle, making this all very, very ugly.
readBamSequence :: FilePath -> BamIndex -> Refseq -> IO [ BamRec ]
readBamSequence fp idx refseq = do
    case idx ! refseq of
        0 -> return []
        virtoff -> do
            let uoffset = virtoff .&. 0xffff
                coffset = virtoff `shiftR` 16
            hdl <- openFile fp ReadMode
            hSeek hdl AbsoluteSeek (fromIntegral coffset)
            go . L.drop (fromIntegral uoffset) . decompressBgzf <$> L.hGetContents hdl
  where
    go s = case runGetState getBamEntry' s 0 of
            _ | L.null s                           -> []
            (Just a, s',  _) | b_rname a == refseq -> a : go s'
                             | otherwise           -> []
            (Nothing, s', _)                       -> go s'

    getBamEntry' = isEmpty >>= \e -> if e then return Nothing else Just <$> getBamEntry
-}



-- | A list of reference sequences.  Note that the first reference sequence must have index 0
type Refs = Array Refseq (Seqid, Int)
type MRefs = IOArray Refseq (Seqid, Int)

-- | The empty list of references.  Needed for BAM files that don't really store alignments.
noRefs :: Refs
noRefs = listArray (toEnum 1, toEnum 0) []



-- | Decodes a raw block into a @BamRec@.
decodeBamEntry :: BamRaw -> BamRec
decodeBamEntry (BamRaw offs s) = case G.runGet go s of
    (Left  e,  _)             -> error e
    (Right r, s') | S.null s' -> r
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
    w2c :: Word8 -> Char
    w2c = chr . fromIntegral

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
-- We start new BGZF blocks at sensible places, then return them
-- individually.  Concatening the result gives a valid file, but the
-- chunks can be turned into an index, too.
encodeBam :: Monad m => BamMeta -> Refs -> Iteratee S.ByteString m a -> Iteratee S.ByteString m a
encodeBam meta refs i = lift (I.enumPure1Chunk header i >>= I.enumPure1Chunk S.empty) >>= compress
  where 
    header = S.concat . L.toChunks $ runPut putHeader

    putHeader = do putWord32le 0x42414D01 -- "BAM\1"
                   let hdr = showBamMeta meta L.empty
                   putWord32le $ fromIntegral $ L.length hdr
                   putLazyByteString hdr
                   put_int_32 . rangeSize $ bounds refs
                   mapM_ putRef $ elems refs

    putRef (n,l) = do put_int_32 $ S.length n + 1
                      putByteString n
                      putWord8 0
                      put_int_32 l


encodeBamEntry :: BamRec -> S.ByteString
encodeBamEntry = S.concat . L.toChunks . runPut . putEntry    
  where
    putEntry e = do let c = runPut $ putEntry' e
                    put_int_32 $ L.length c
                    putLazyByteString c

    putEntry' b = do putWord32le   $ unRefseq $ b_rname b
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
                     forM_ (M.toList $ b_exts b) $ \(k,v) -> 
                        case k of [c,d] -> putChr c >> putChr d >> putValue v
                                  _     -> error $ "invalid field key " ++ show k

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
-- file.  We also keep the offsets of those blocks and write a separate
-- binary index (for Sector).
-- XXX Do we need this?  Do we want it?
-- Note: The index format of Sector is underspecified.  We write the
-- offset as little endian, but it's probably platform dependent.
{-
writeBamFile :: FilePath -> Bam -> IO ()
writeBamFile fp bam =
    withBinaryFile fp WriteMode                 $ \h_bam ->
    withBinaryFile (fp ++ ".idx") WriteMode     $ \h_idx ->
    let go _ [    ] = return ()                                     -- done
        go p (b:bs) = do
            L.hPut h_bam b                                          -- write block
            L.hPut h_idx $ runPut $ putWord64le (fromIntegral p)    -- to index, too
            go (p+L.length b) bs                                    -- recurse
    in case encodeBam bam of
        []     -> return ()                 -- can't actually happen
        (b:bs) -> do L.hPut h_bam b         -- write header, not part of index
                     go (L.length b) bs     -- recurse
-}

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
 flagAuxillary, flagFailsQC, flagDuplicate, flagLowQuality, flagLowComplexity :: Int

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
flagLowQuality = 0x800
flagLowComplexity = 0x1000


isPaired, isProperlyPaired, isUnmapped, isMateUnmapped, isReversed, isMateReversed, isAuxillary, isFailsQC, isDuplicate,
 isLowQuality, isLowComplexity :: BamRec -> Bool

isPaired         = flip testBit  0 . b_flag
isProperlyPaired = flip testBit  1 . b_flag
isUnmapped       = flip testBit  2 . b_flag
isMateUnmapped   = flip testBit  3 . b_flag
isReversed       = flip testBit  4 . b_flag
isMateReversed   = flip testBit  5 . b_flag
isAuxillary      = flip testBit  8 . b_flag
isFailsQC        = flip testBit  9 . b_flag
isDuplicate      = flip testBit 10 . b_flag
isLowQuality     = flip testBit 11 . b_flag
isLowComplexity  = flip testBit 12 . b_flag

 
-- | tests if a record is a "first mate"
-- Returns true if the read is flagged as "paired" and "first mate".
-- Does not return true for single reads.
isFirstMate :: BamRec -> Bool
isFirstMate = (==) good . (.&.) type_mask . b_flag
  where good = flagPaired .|. flagFirstMate

-- | Tests if a record is a "second mate".
-- Returns true if the read is flagged as "paired" and "second mate".
-- Does not return true for single reads, even if they are adapter
-- trimmed or merged.
isSecondMate :: BamRec -> Bool
isSecondMate = (==) good . (.&.) type_mask . b_flag
  where good = flagPaired .|. flagSecondMate

-- | Tests if a read is merged.
-- Returns true for merged single reads (not flagged as "paired", but
-- flagged as both "first mate" and "second mate"), does not get
-- confused by true singles (not flagged as "first mate" and "second
-- mate").
isMerged :: BamRec -> Bool
isMerged = (==) good . (.&.) type_mask . b_flag
  where good = flagFirstMate .|. flagSecondMate

-- | Tests if a read is adapter trimmed.
-- Returns true for adapter trimmed single reads (not flagged as
-- "paired", but flagged as "second mate"), does not get confused by
-- paired reads (also flagged as "paired") or merged reads (also flagged
-- as "first mate").
isAdapterTrimmed :: BamRec -> Bool
isAdapterTrimmed = (==) flagSecondMate . (.&.) type_mask . b_flag


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

compareNames :: L.ByteString -> L.ByteString -> Ordering
compareNames n m = case (L.uncons n, L.uncons m) of
        ( Nothing, Nothing ) -> EQ
        ( Just  _, Nothing ) -> GT
        ( Nothing, Just  _ ) -> LT
        ( Just (c,n'), Just (d,m') )
            | isDigit c && isDigit d -> 
                let Just (u,n'') = L.readInt n
                    Just (v,m'') = L.readInt m
                in case u `compare` v of 
                    LT -> LT
                    GT -> GT
                    EQ -> n'' `compareNames` m''
            | otherwise -> case c `compare` d of 
                    LT -> LT
                    GT -> GT
                    EQ -> n' `compareNames` m'
                                         

-- | Gets all the filter flags.
-- Filter flags are the standard "fails QC", "is duplicate" and our
-- extensions "fails quality" and "fails complexity".
filterFlags :: BamRec -> Int
filterFlags = (.&.) mask . b_flag  
  where mask = flagLowQuality .|. flagLowComplexity .|. flagFailsQC .|. flagDuplicate

-- Stop gap solution.  we only get the first offset from the linear
-- index, which allows us to navigate to a target sequence.  Will do the
-- rest when I need it.
type BamIndex = Array Refseq Word64

readBamIndex :: L.ByteString -> BamIndex
readBamIndex = runGet getBamIndex
  where
    getBamIndex = do
        magic <- getByteString 4
        when (magic /= "BAI\1" ) $ fail "BAI signature not found"
        nref <- get_int_32
        offs <- replicateM nref $ do
                    nbins <- get_int_32
                    replicateM_ nbins $ do
                        _bin <- get_int_32 -- "distinct bin", whatever that means
                        nchunks <- get_int_32
                        replicateM_ nchunks $ getWord64le >> getWord64le
                    nintv <- get_int_32
                    os <- filter (/= 0) <$> replicateM nintv getWord64le
                    return $ if null os then 0 else minimum os
        return $! listArray (toEnum 0, toEnum (nref-1)) offs

    get_int_32 = fromIntegral <$> getWord32le


data BamMeta = BamMeta {
        meta_hdr :: BamHeader,
        meta_other_shit :: [(Char, Char, BamOtherShit)],
        meta_comment :: [S.ByteString]
    } deriving Show

data BamHeader = BamHeader {
        hdr_version :: (Int, Int),
        hdr_sorting :: BamSorting,
        hdr_other_shit :: BamOtherShit
    } deriving Show

data BamSorting = Unsorted | Grouped | Queryname | Coordinate | GroupSorted 
    deriving Show

type BamOtherShit = [(Char, Char, S.ByteString)]

parseBamMeta :: P.Parser BamMeta
parseBamMeta = foldr ($) nullMeta <$> P.many pLine 
  where
    pLine = P.char '@' >> P.choice [hdLine, coLine, otherLine] <* P.char '\n'
    hdLine = P.string "HD\t" >> 
             (\fns meta -> meta { meta_hdr = foldr ($) (meta_hdr meta) fns })
               <$> P.sepBy1 (P.choice [hdvn, hdso, hdother]) (P.char '\t')
    
    hdvn = P.string "VN:" >>
           (\a b hdr -> hdr { hdr_version = (a,b) })
             <$> P.decimal <*> (P.char '.' >> P.decimal)

    hdso = P.string "SO:" >>
           (\s hdr -> hdr { hdr_sorting = s })
             <$> P.choice [ Grouped  <$ P.string "grouped"
                          , Queryname <$ P.string "queryname"
                          , Coordinate <$ P.string "coordinate"
                          , GroupSorted <$ P.string "groupsort"
                          , Unsorted <$ P.skipWhile (\c -> c/='\t' && c/='\n') ]

    hdother = (\t hdr -> hdr { hdr_other_shit = t : hdr_other_shit hdr }) <$> tagother
    
    coLine = P.string "CO\t" >>
             (\s meta -> meta { meta_comment = s : meta_comment meta })
               <$> P.takeWhile (/= 'n')

    otherLine = (\a b ts meta -> meta { meta_other_shit = (a,b,ts) : meta_other_shit meta })
                  <$> P.anyChar <*> P.anyChar <*> (P.char '\t' >> P.sepBy1 tagother (P.char '\t'))

    tagother :: P.Parser (Char,Char,S.ByteString)
    tagother = (,,) <$> P.anyChar <*> P.anyChar <*> (P.char ':' >> P.takeWhile (\c -> c/='\t' && c/='\n'))

nullMeta :: BamMeta
nullMeta = BamMeta (BamHeader (0,0) Unsorted []) [] []

showBamMeta :: BamMeta -> L.ByteString -> L.ByteString
showBamMeta (BamMeta h o c) = 
    show_bam_meta_hdr h .
    foldr ((.) . show_bam_meta_other) id o .
    foldr ((.) . show_bam_meta_comment) id c
  where
    show_bam_meta_hdr (BamHeader (major,minor) so os) = 
        L.append "@HD\tVN:" . L.append (L.pack (show major ++ ':' : show minor)) .
        L.append (case so of Unsorted -> L.empty
                             Grouped  -> "\tSO:grouped"
                             Queryname  -> "\tSO:queryname"
                             Coordinate  -> "\tSO:coordinate"
                             GroupSorted  -> "\tSO:groupsort") .
        foldr ((.) . show_bam_other) id os .
        L.cons '\n'

    show_bam_meta_comment cm = L.append "@CO\t" . L.append (L.fromChunks [cm]) . L.cons '\n'

    show_bam_meta_other (a,b,ts) = 
        L.cons '@' . L.cons a . L.cons b . 
        foldr ((.) . show_bam_other) id ts . L.cons '\n'

    show_bam_other (a,b,v) = 
        L.cons '\t' . L.cons a . L.cons b . L.cons ':' . L.append (L.fromChunks [v])


-- | Bam record in its native encoding along with virtual address.
data BamRaw = BamRaw { virt_offset :: {-# UNPACK #-} !Int64
                     , raw_data :: !S.ByteString }


-- | Decode a BAM stream into raw entries.  Note that the entries can be
-- unpacked using @decodeBamEntry@.
decodeBam :: MonadIO m => (BamMeta -> Refs -> Iteratee [BamRaw] m a) -> Iteratee Block m a
decodeBam inner = do meta <- liftBlock get_bam_header
                     refs <- liftBlock get_ref_array
                     loop $ inner meta refs
  where
    loop it = I.isFinished >>= loop' it
    loop' it True = lift $ run it
    loop' it False = do off <- getOffset
                        it' <- liftBlock $ do bsize <- endianRead4 LSB
                                              raw <- getString (fromIntegral bsize)
                                              lift $ I.enumPure1Chunk [BamRaw off raw] it
                        loop it'

    get_bam_header  = do magic <- I.heads "BAM\SOH"
                         when (magic /= 4) $ fail "BAM signature not found"
                         hdr_len <- endianRead4 LSB
                         I.joinI $ I.take (fromIntegral hdr_len) $ parserToIteratee parseBamMeta

    get_ref_array = do nref <- endianRead4 LSB
                       let lr = (toEnum 0,toEnum (fromIntegral nref-1))
                       arr <- liftIO $ newArray_ lr
                       forM_ (range lr) $ \i -> do namelength <- endianRead4 LSB
                                                   nm <- getString (fromIntegral namelength)
                                                   ln <- endianRead4 LSB
                                                   liftIO $ writeArray arr i (S.init nm,fromIntegral ln)
                       liftIO $ unsafeFreeze (arr :: MRefs)


-- ------------------------------------------------------------------- Tests

test :: IO ()
test = I.fileDriver (decompress' (decodeBam (\_ _ -> print_names))) 
                    "/mnt/ngs_data/101203_SOLEXA-GA04_00007_PEDi_MM_QF_SR/Ibis/BWA/s_5_L3280_sequence_mq_hg19_nohap.bam"

print_names :: Iteratee [BamRaw] IO ()
print_names = I.mapM_ $ S.putStrLn . b_qname . decodeBamEntry 
    
