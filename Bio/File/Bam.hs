module Bio.File.Bam (
    module Bio.Base,

    decompressBgzf,
    compressBgzfDynamic,
    compressBgzfSingle,
    isBam,
    readBam,
    encodeBam,
    writeBamFile,

    CodedSeq(..),
    CodedCigar(..),
    Bam(..),
    BamRec(..),
    Refs,
    noRefs,
    Ext(..),
    Refseq(..),
    invalidRefseq,
    invalidPos,
    isValidRefseq,
    compareNames,

    MdOp(..),
    getMd,
    readMd,

    CigOp(..),
    cigOp,
    cigLen,
    mkCigOp,
    cigarToAlnLen,
    mkExtKey,

    decodeSeq,
    encodeSeq,
    inflateSeq,
    deflateSeq,
    packCigar,
    unpackCigar,

    bamFlagPaired,
    bamFlagProperlyPaired,
    bamFlagUnmapped,
    bamFlagMateUnmapped,
    bamFlagStrand,
    bamFlagMateStrand,
    bamFlagFirstMate,
    bamFlagSecondMate,
    bamFlagAuxillary,
    bamFlagFailsQC,
    bamFlagDuplicate,
    bamFlagQualityFiltered,
    bamFlagComplexityFiltered,

    bamFilterFlags,
    isFirstMate,
    isSecondMate,
    isMerged,
    isAdapterTrimmed
) where

import Bio.Base

import Codec.Compression.GZip
import Control.Monad ( replicateM )
import Control.Applicative ((<$>), (<*>), (*>) )
import Data.Array.Unboxed
import Data.Binary.Get
import Data.Binary.Put
import Data.Bits ( testBit, shiftL, shiftR, (.&.), (.|.) )
import Data.Char ( chr, ord, isDigit )
import Data.Int  ( Int64 )
import Data.List ( genericTake, genericLength )
import Data.Word ( Word32, Word8 )
import System.IO ( withBinaryFile, IOMode(WriteMode) )

import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.ByteString.Lazy as LB
import qualified Data.ByteString.Char8 as S
import qualified Data.Map as M

-- | Sequence in BAM coding
-- Bam encodes two bases per byte, we store the length separately to
-- allow reconstruction without trailing junk.
data CodedSeq = CodedSeq { coded_seq_length :: !Int64, coded_seq_bases :: L.ByteString }

instance Show CodedSeq where show = show . decodeSeq

-- | Cigar line in BAM coding
-- Bam encodes an operation and a length into a single integer, we keep
-- those integers in an array.
data CodedCigar = CodedCigar { unCodedCigar :: UArray Int Int }

instance Show CodedCigar where show = show . elems . unCodedCigar

-- | extracts the aligned length from a cigar line
-- This gives the length of an alignment as measured on the reference,
-- which is different from the length on the query or the length of the
-- alignment.
cigarToAlnLen :: CodedCigar -> Int
cigarToAlnLen (CodedCigar cig) = sum $ map l $ elems cig
  where l c = let op = c .&. 0xf
              in if op == 0 || op == 2 || op == 3 then c `shiftR` 4 else 0
    

-- | Reference sequence in Bam
-- Bam enumerates the reference sequences and then sorts by index.  We
-- need to track that index if we want to reproduce the sorting order.
newtype Refseq = Refseq { unRefseq :: Word32 } deriving (Show, Eq, Ord)

-- | tests whether a reference sequence is valid
-- Returns true unless the the argument equals 'invalidRefseq'.
isValidRefseq :: Refseq -> Bool
isValidRefseq = (/=) invalidRefseq

-- | the invalid Refseq
-- Bam uses this value to encode an invalid reference sequence.
invalidRefseq :: Refseq
invalidRefseq = Refseq 0xffffffff

invalidPos :: Int
invalidPos = -1

-- | internal representation of a BAM record
data BamRec = BamRec {
    b_qname :: !L.ByteString,
    b_flag  :: !Int,
    b_rname :: !Refseq,
    b_pos   :: !Int,
    b_mapq  :: !Int,
    b_cigar :: !CodedCigar,
    b_mrnm  :: !Refseq,
    b_mpos  :: !Int,
    b_isize :: !Int,
    b_seq   :: !CodedSeq,
    b_qual  :: !L.ByteString,       -- ^ quality, may be empty
    b_exts  :: M.Map Int Ext,
    b_virtual_offset :: !Int64      -- ^ virtual offset for indexing purposes
} deriving Show


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

-- | decompresses Bgzf or Gzip
-- This checks for the Bgzf header, and if present, decompresses Bgzf
-- chunks.  Else it decompresses GZip.
decompressBgzf :: L.ByteString -> L.ByteString
decompressBgzf = go
  where
    go s | L.null s = L.empty 
    go s = case runGet get_bgzf_hdr s of 
                Nothing -> if isGzip s then decompress s else s
                Just l  -> case L.splitAt (fromIntegral l + 1) s of 
                                (u,v) -> decompress u `L.append` go v
                     
    get_bgzf_hdr = do id1 <- getWord8
                      id2 <- getWord8
                      skip 1
                      flg <- getWord8
                      if id1 == 31 && id2 == 139 && flg `testBit` 2 
                        then do skip 6
                                xdata <- getWord16le >>= getLazyByteString . fromIntegral
                                return $ runGet get_bsize xdata
                        else return Nothing

    get_bsize = isEmpty >>= \e -> if e then return Nothing else do
                    i1 <- getWord8
                    i2 <- getWord8
                    len <- getWord16le
                    if i1 == 66 && i2 == 67 && len == 2 
                      then Just <$> getWord16le
                      else skip (fromIntegral len) >> get_bsize

maxBlockSize :: Int64
maxBlockSize = 65450 -- 64k with some room for headers and uncompressible stuff

-- | The EOF marker for BGZF files.
-- This is just an empty string compressed as BGZF.  Append to BAM files
-- to indicate their end.
bgzfEofMarker :: L.ByteString
bgzfEofMarker = compressBgzfSingle L.empty

-- | Compresses a list of strings into a dynamic BGZF stream.
-- We try to create large blocks (no more than 64k, of course), but if
-- we have to start a new block, we try to do it at the start of a new
-- input element.  Inputs are compressed using 'compressBgzfSingle',
-- possibly resulting in new blocks.  The output can be concatenated
-- into a BGZF stream, output elements correspond to block boundaries
-- that coincide with input block boundaries.
compressBgzfDynamic :: [ L.ByteString ] -> [ L.ByteString ]
compressBgzfDynamic = go 0 []
  where go len acc [    ]                                    = [ compressBgzfSingle $ L.concat $ reverse acc ]
        go len acc (s:ss) | len + L.length s <= maxBlockSize = go (len + L.length s) (s:acc) ss
                          | null acc                         = go (L.length s) [s] ss
                          | otherwise                        = (compressBgzfSingle $ L.concat $ reverse acc) : go (L.length s) [s] ss
    
-- | Compress a single string into a BGZF stream.
-- The resulting stream has arbitrary blocks, as the maximum size
-- dictates, and no EOF marker.  Technically the result is a valid file,
-- but this more valuable as a building block.
compressBgzfSingle :: L.ByteString -> L.ByteString
compressBgzfSingle s = hdr `L.append` rest `L.append` cont
  where
    (l,r) = L.splitAt maxBlockSize s
    cont  = if L.null r then L.empty else compressBgzfSingle r
    z = compress l
    (hdr, rest, _) = runGetState patch_header z 0
    patch_header = do k <- getWord16le
                      m <- getWord8
                      f <- getWord8
                      t <- getWord32le
                      xf <- getWord8
                      getWord8 -- OS
                      xlen <- if f `testBit` 2 then getWord16le else return 0

                      return $ runPut $ do 
                            putWord16le k
                            putWord8 m
                            putWord8 $ f .|. 4
                            putWord32le t
                            putWord8 xf
                            putWord8 0xff
                            putWord16le $ xlen + 6
                            putWord8 66
                            putWord8 67
                            putWord16le 2
                            putWord16le . fromIntegral $ L.length z + 5 + 
                                if f `testBit` 2 then 0 else 2


readBam :: L.ByteString -> Bam
readBam s0 = Bam hdr refs (go s1 p1)
  where
    ((hdr,refs), s1, p1) = runGetState getBamHeader (decompressBgzf s0) 0

    go s p = case runGetState getBamEntry' s p of
            _ | L.null s      -> []
            (Just a, s',  p') -> a : go s' p'
            (Nothing, s', p') -> go s' p'

    getBamEntry' = isEmpty >>= \e -> if e then return Nothing else Just <$> getBamEntry

-- note: first reference sequence must have index 0
type Refs = Array Int (Seqid, Int)

data Bam = Bam { bam_header :: L.ByteString
               , bam_refs   :: Refs
               , bam_recs   :: [ BamRec ] } deriving Show

noRefs :: Refs
noRefs = listArray (1,0) []

getBamHeader :: Get (L.ByteString,Refs)
getBamHeader = do "BAM\SOH" <- L.unpack <$> getLazyByteString 4
                  hdr <- get_int_32 >>= getLazyByteString
                  nref <- get_int_32
                  let loop acc 0 = return $! array (0,nref-1) $ zip [nref-1, nref-2 ..] acc
                      loop acc n = do nm <- get_int_32 >>= fmap S.init . getByteString
                                      ln <- get_int_32
                                      loop ((nm,ln):acc) $! n-1
                  (,) hdr <$> loop [] nref                                  


get_int_32, get_int_16, get_int_8 :: Integral a => Get a
get_int_32 = fromIntegral <$> getWord32le
get_int_16 = fromIntegral <$> getWord16le
get_int_8  = fromIntegral <$> getWord8


getBamEntry :: Get BamRec
getBamEntry = do
    start_offs <- bytesRead
    end_pos <- (+) <$> (get_int_32 :: Get Int64) <*> bytesRead 
    rid <- Refseq <$> getWord32le
    start <- get_int_32
    skip 1  -- name length
    mapq <- get_int_8
    skip 2  -- bin number
    cigar_len <- get_int_16
    flag <- get_int_16
    read_len <- get_int_32
    mate_rid <- Refseq <$> getWord32le
    mate_pos <- get_int_32
    ins_size <- get_int_32

    read_name <- getLazyByteStringNul
    cigar <- listArray (1,cigar_len) <$> replicateM cigar_len get_int_32
    qry_seq <- getLazyByteString ((read_len+1) `div` 2)
    qual <- (\qs -> if LB.all (0xff ==) qs then L.empty else qs) <$> getLazyByteString read_len
    exts <- getExtensions end_pos M.empty

    return $ BamRec read_name flag rid start mapq (CodedCigar cigar)
                    mate_rid mate_pos ins_size (CodedSeq read_len qry_seq) qual exts start_offs
        

data Ext = Int Int | Text L.ByteString | Char Word8 deriving Show

getExtensions :: Int64 -> M.Map Int Ext -> Get (M.Map Int Ext)
getExtensions p m = do p' <- bytesRead
                       case p `compare` p' of
                            GT -> getExt >>= \(k,a) -> getExtensions p (M.insert k a m)
                            EQ -> return m
                            LT -> error "corrupt file or bug"

getExt :: Get (Int, Ext)
getExt = do key <- get_int_16
            typ <- getWord8
            res <- case chr . fromIntegral $ typ of
                    'c' -> Int <$> get_int_8
                    'C' -> Int <$> get_int_8
                    's' -> Int <$> get_int_16
                    'S' -> Int <$> get_int_16
                    'i' -> Int <$> get_int_32
                    'I' -> Int <$> get_int_32
                    'Z' -> Text <$> getLazyByteStringNul
                    'A' -> Char <$> getWord8
                    x   -> error $ "cannot handle optional field type " ++ [x]
            return (key,res)

mkExtKey :: Char -> Char -> Int
mkExtKey x y = ord x .|. (ord y `shiftL` 8)

-- | Encode stuff into a BAM stream.
-- We start new BGZF blocks at sensible places, then return them
-- individually.  Concatening the result gives a valid file, but the
-- chunks can be turned into an index, too.
encodeBam :: Refs -> [ BamRec ] -> [ L.ByteString ]
encodeBam refs xs = compressBgzfSingle (runPut putHeader) :
                    compressBgzfDynamic (map (runPut . putEntry) xs) ++ 
                    [ bgzfEofMarker ]
  where 
    putHeader = do putByteString $ S.pack "BAM\SOH"
                   putWord32le 0
                   let (l,r) = bounds refs 
                   put_int_32 (r-l+1)
                   mapM_ putRef $ elems refs

    putRef (n,l) = do put_int_32 $ S.length n + 1
                      putByteString n
                      putWord8 0
                      put_int_32 l

    putEntry e = do let c = runPut $ putEntry' e
                    put_int_32 $ L.length c
                    putLazyByteString c

    putEntry' b = do putWord32le $ unRefseq $ b_rname b
                     put_int_32 $ b_pos b
                     put_int_8 $ L.length (b_qname b) + 1
                     put_int_8 $ b_mapq b
                     put_int_16 $ bin (b_pos b) (b_pos b + cigarToAlnLen (b_cigar b) - 1)
                     let (l,r) = bounds $ unCodedCigar $ b_cigar b
                     put_int_16 (r-l+1)
                     put_int_16 $ b_flag b
                     put_int_32 $ coded_seq_length $ b_seq b
                     putWord32le $ unRefseq $ b_mrnm b
                     put_int_32 $ b_mpos b
                     put_int_32 $ b_isize b
                     putLazyByteString $ b_qname b
                     putWord8 0
                     mapM_ put_int_32 $ elems $ unCodedCigar $ b_cigar b
                     putLazyByteString $ coded_seq_bases $ b_seq b
                     putLazyByteString $ if not (L.null (b_qual b)) then b_qual b
                                         else LB.replicate (coded_seq_length $ b_seq b) 0xff
                     mapM_ (\(k,v) -> put_int_16 k >> putValue v) $ M.toList $ b_exts b
                     
-- | writes stuff to a BAM file                     
-- We generate BAM with dynamic blocks, then stream them out to the
-- file.  We also keep the offsets of those blocks and write a separate
-- binary index (for Sector).
-- Note: The index format of Sector is underspecified.  We write the
-- offset as little endian, but it's probably platform dependent.
writeBamFile :: FilePath -> Refs -> [ BamRec ] -> IO ()
writeBamFile fp refs alns =
    withBinaryFile fp WriteMode                 $ \h_bam ->
    withBinaryFile (fp ++ ".idx") WriteMode     $ \h_idx ->
    let go _ [    ] = return ()                                     -- done
        go p (b:bs) = do
            L.hPut h_bam b                                          -- write block
            L.hPut h_idx $ runPut $ putWord64le (fromIntegral p)    -- to index, too
            go (p+L.length b) bs                                    -- recurse
    in case encodeBam refs alns of
        []     -> return ()                 -- can't actually happen
        (b:bs) -> do L.hPut h_bam b         -- write header, not part of index
                     go (L.length b) bs     -- recurse


put_int_32, put_int_16, put_int_8 :: Integral a => a -> Put
put_int_32 = putWord32le . fromIntegral
put_int_16 = putWord16le . fromIntegral
put_int_8  = putWord8 . fromIntegral

putChr :: Char -> Put
putChr = putWord8 . fromIntegral . ord

bin :: Int -> Int -> Int
bin beg end = mkbin 14 $ mkbin 17 $ mkbin 20 $ mkbin 23 $ mkbin 16 $ 0
  where mkbin n x = if beg `shiftR` n /= end `shiftR` n then x
                    else ((1 `shiftL` (29-n))-1) `div` 7 + (beg `shiftR` n)


putValue :: Ext -> Put
putValue (Text t) = putChr 'Z' >> putLazyByteString t >> putWord8 0
putValue (Char c) = putChr 'A' >> putWord8 c
putValue (Int i) | i < -0xffff = putChr 'i' >> put_int_32 i
                 | i < -0xff   = putChr 's' >> put_int_16 i
                 | i < 0       = putChr 'c' >> put_int_8  i
                 | i > -0xffff = putChr 'I' >> put_int_32 i
                 | i > -0xff   = putChr 'S' >> put_int_16 i
                 | otherwise   = putChr 'C' >> put_int_8  i

data CigOp = Mat | Ins | Del | Nop | SMa | HMa | Pad 
    deriving ( Eq, Ord, Enum, Show, Bounded, Ix )

cigOp :: Int -> CigOp
cigOp c | cc <= fromEnum (maxBound :: CigOp) = toEnum cc
        | otherwise = error "unknown Cigar operation"
  where cc = c .&. 0xf

cigLen :: Int -> Int
cigLen c = c `shiftR` 4

mkCigOp :: CigOp -> Int -> Int
mkCigOp op len = (fromIntegral len `shiftL` 4) .|. fromEnum op

-- | unpack a two-bases-per-byte sequence
inflateSeq :: CodedSeq -> [Word8]
inflateSeq (CodedSeq l s) = genericTake l $ go s
  where
    go s | L.null s = []
         | otherwise = let x = LB.head s in x `shiftR` 4 : x .&. 0xf : go (L.tail s)

-- | repeatedly packs two bases into one byte
deflateSeq :: [Word8] -> CodedSeq
deflateSeq ws = CodedSeq (genericLength ws) (LB.pack $ go ws)
  where
    go (x:y:zs) = ((x `shiftL` 4) .|. y) : go zs
    go [x] = [x `shiftL` 4]
    go [] = []

-- | converts a Bam sequence into Nucleotides
-- Everything that isn't representable maps to N.
decodeSeq :: CodedSeq -> [Nucleotide]
decodeSeq = map (bases !) . inflateSeq
  where bases = listArray (0,15) (map toNucleotide "NACNGNNNTNNNNNNN") :: Array Word8 Nucleotide

-- | converts a sequence of Nucleotides into the Bam representaion
encodeSeq :: [Nucleotide] -> CodedSeq
encodeSeq = deflateSeq . map num
  where
    num A = 1
    num C = 2
    num G = 4
    num T = 8
    num N = 15
    num _ = 0



packCigar :: [Int] -> CodedCigar
packCigar cs = CodedCigar $ listArray (1, length cs) cs

unpackCigar :: CodedCigar -> [Int]
unpackCigar = elems . unCodedCigar


data MdOp = MdNum Int | MdRep Nucleotide | MdDel [Nucleotide] deriving Show

getMd :: BamRec -> Maybe [MdOp]
getMd r = case M.lookup (mkExtKey 'M' 'D') $ b_exts r of
    Just (Text mdfield) -> readMd mdfield
    _                   -> Nothing

readMd :: L.ByteString -> Maybe [MdOp]
readMd s | L.null s           = return []
         | isDigit (L.head s) = do (n,t) <- L.readInt s
                                   (MdNum n :) <$> readMd t
         | L.head s == '^'    = let (a,b) = L.break isDigit (L.tail s)
                                in (MdDel (map toNucleotide $ L.unpack a) :) <$> readMd b
         | otherwise          = (MdRep (toNucleotide $ L.head s) :) <$> readMd (L.tail s)


bamFlagPaired, bamFlagProperlyPaired, bamFlagUnmapped,
 bamFlagMateUnmapped, bamFlagStrand, bamFlagMateStrand,
 bamFlagFirstMate, bamFlagSecondMate, bamFlagAuxillary,
 bamFlagFailsQC, bamFlagDuplicate :: Int

bamFlagPaired = 0x1
bamFlagProperlyPaired = 0x2
bamFlagUnmapped = 0x4
bamFlagMateUnmapped = 0x8
bamFlagStrand = 0x10
bamFlagMateStrand = 0x20
bamFlagFirstMate = 0x40
bamFlagSecondMate = 0x80
bamFlagAuxillary = 0x100
bamFlagFailsQC = 0x200
bamFlagDuplicate = 0x400

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
                                         

-- * Bam functions depending on MPI conventions

bamFlagQualityFiltered, bamFlagComplexityFiltered :: Int
bamFlagQualityFiltered = 0x800
bamFlagComplexityFiltered = 0x1000

-- | get all the filter flags
-- Filter flags are the standard "fails QC", "is duplicate" and our
-- extensions "fails quality" and "fails complexity".
bamFilterFlags :: BamRec -> Int
bamFilterFlags = (.&. mask) . b_flag  
  where mask = bamFlagQualityFiltered .|. bamFlagComplexityFiltered .|. bamFlagFailsQC

-- | tests if a record is a "first mate"
-- Returns true if the read is flagged as "paired" and "first mate".
-- Does not return true for single reads.
isFirstMate :: BamRec -> Bool
isFirstMate = (== good) . (.&. mask) . b_flag
  where good = bamFlagPaired .|. bamFlagFirstMate
        mask = good .|. bamFlagSecondMate

-- | tests if a record is a "second mate"
-- Returns true if the read is flagged as "paired" and "second mate".
-- Does not return true for single reads, even if they are adapter
-- trimmed.
isSecondMate :: BamRec -> Bool
isSecondMate = (== good) . (.&. mask) . b_flag
  where good = bamFlagPaired .|. bamFlagSecondMate
        mask = good .|. bamFlagFirstMate

-- | tests if a read is merged
-- Returns true for merged single reads, does not get confused by true
-- singles.
isMerged :: BamRec -> Bool
isMerged = (== good) . (.&. mask) . b_flag
  where good = bamFlagFirstMate .|. bamFlagSecondMate
        mask = good

-- | tests if a read is adapter trimmed
-- Returns true for adapter trimmed single reads, does not get confused
-- by paired reads.
isAdapterTrimmed :: BamRec -> Bool
isAdapterTrimmed = (== good) . (.&. mask) . b_flag
  where good = bamFlagSecondMate
        mask = bamFlagFirstMate .|. bamFlagPaired .|. good

