module Bio.File.Bam (
    decompress_bgzf,
    compress_bgzf,
    is_bam,
    readBam,
    writeBam,

    CodedSeq(..),
    CodedCigar(..),
    BamRec(..),
    Refs,
    Ext(..),
    Refseq(..),
    invalidRefseq,
    isValidRefseq,

    MdOp(..),
    getMd,

    cig_op,
    cig_len,
    mk_cig_op,
    cigar_to_aln_len,
    mk_ext_key,

    decode_seq,
    inflate_seq,
    pack_cigar

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
import Data.Int ( Int64 )
import Data.List ( genericTake )
import Data.Word ( Word32, Word8 )
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.ByteString.Lazy as LB
import qualified Data.ByteString.Char8 as S
import qualified Data.Map as M

-- | Sequence in BAM coding
-- Bam encodes two bases per byte, we store the length separately to
-- allow reconstruction without trailing junk.
data CodedSeq = CodedSeq { coded_seq_length :: !Int64, coded_seq_bases :: L.ByteString }

instance Show CodedSeq where show = show . decode_seq

-- | Cigar line in BAM coding
-- Bam encodes an operation and a length into a single integer, we keep
-- those integers in an array.
data CodedCigar = CodedCigar { unCodedCigar :: UArray Int Int }

instance Show CodedCigar where show = show . elems . unCodedCigar

-- | extracts the aligned length from a cigar line
-- This gives the length of an alignment as measured on the reference,
-- which is different from the length on the query or the length of the
-- alignment.
cigar_to_aln_len :: CodedCigar -> Int
cigar_to_aln_len (CodedCigar cig) = sum $ map l $ elems cig
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
    b_qual  :: !L.ByteString,
    b_exts  :: M.Map Int Ext,
    b_virtual_offset :: !Int64
} deriving Show


is_bam :: L.ByteString -> Bool
is_bam s = L.pack "BAM\SOH" `L.isPrefixOf` decompress_bgzf s

decompress_bgzf :: L.ByteString -> L.ByteString
decompress_bgzf = go
  where
    go s = case runGet get_bgzf_hdr s of 
                _ | L.null s -> L.empty 
                Nothing -> s
                Just l  -> case L.splitAt (fromIntegral l + 1) s of 
                    (u,v) -> decompress u `L.append` go v
                     
    get_bgzf_hdr = do id1 <- getWord8
                      id2 <- getWord8
                      skip 1
                      flg <- getWord8
                      skip 6
                      xdata <- getWord16le >>= getLazyByteString . fromIntegral
                      if id1 == 31 && id2 == 139 && flg `testBit` 2 
                        then return $ runGet get_bsize xdata
                        else return Nothing

    get_bsize = isEmpty >>= \e -> if e then return Nothing else do
                    i1 <- getWord8
                    i2 <- getWord8
                    len <- getWord16le
                    if i1 == 66 && i2 == 67 && len == 2 
                      then Just <$> getWord16le
                      else skip (fromIntegral len) >> get_bsize

compress_bgzf :: L.ByteString -> L.ByteString
compress_bgzf s | L.null s = s 
                | otherwise = hdr `L.append` rest `L.append` compress_bgzf r
  where
    (l,r) = L.splitAt 65000 s
    z = compress l
    (hdr, rest, _) = runGetState patch_header z 0
    patch_header = do k <- getWord16le
                      m <- getWord8
                      f <- getWord8
                      t <- getWord32le
                      xf <- getWord16le
                      xlen <- if f `testBit` 2 then getWord16le else return 0

                      return $ runPut $ do 
                            putWord16le k
                            putWord8 m
                            putWord8 $ f .|. 4
                            putWord32le t
                            putWord16le xf
                            putWord16le $ xlen + 6
                            putWord8 66
                            putWord8 67
                            putWord16le 2
                            putWord16le . fromIntegral $ L.length z + 5 + 
                                if f `testBit` 2 then 0 else 2


readBam :: L.ByteString -> ( Refs, [ BamRec ] )
readBam s0 = ( refs, go s1 p1 )
  where
    (refs, s1, p1) = runGetState getBamHeader (decompress_bgzf s0) 0

    go s p = case runGetState getBamEntry' s p of
            _ | L.null s      -> []
            (Just a, s',  p') -> a : go s' p'
            (Nothing, s', p') -> go s' p'

    getBamEntry' = isEmpty >>= \e -> if e then return Nothing else Just <$> getBamEntry

type Refs = Array Int (Seqid, Int)

getBamHeader :: Get Refs
getBamHeader = do "BAM\SOH" <- L.unpack <$> getLazyByteString 4
                  get_int_32 >>= skip
                  nref <- get_int_32
                  listArray (0,nref-1) <$> 
                    replicateM nref ((,) <$> (get_int_32 >>= (fmap S.init . getByteString)) <*> get_int_32)

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
    qual <- getLazyByteString read_len
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

mk_ext_key :: Char -> Char -> Int
mk_ext_key x y = ord x .|. (ord y `shiftL` 8)

writeBam :: Refs -> [ BamRec ] -> L.ByteString
writeBam refs xs = compress_bgzf $ runPut $ putHeader >> mapM_ putEntry xs
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
                     put_int_16 $ bin (b_pos b) (b_pos b + cigar_to_aln_len (b_cigar b) - 1)
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
                     putLazyByteString $ b_qual b
                     mapM_ (\(k,v) -> put_int_16 k >> putValue v) $ M.toList $ b_exts b
                     
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

cig_op :: Int -> Int
cig_op c = c .&. 0xf

cig_len :: Int -> Int
cig_len c = c `shiftR` 4

mk_cig_op :: Int -> Int -> Int
mk_cig_op op len = (fromIntegral len `shiftL` 4) .|. op

inflate_seq :: L.ByteString -> [Word8]
inflate_seq s | L.null s = []
              | otherwise = let x = LB.head s in x `shiftR` 4 : x .&. 0xf : inflate_seq (L.tail s)

decode_seq :: CodedSeq -> String
decode_seq (CodedSeq l s) = map (bases !) $ genericTake l $ inflate_seq s
  where bases = listArray (0,15) "=AC.G...T......N" :: UArray Word8 Char

pack_cigar :: [Int] -> CodedCigar
pack_cigar cs = CodedCigar $ listArray (1, length cs) cs


data MdOp = MdNum Int | MdRep Char | MdDel [Char] deriving Show

getMd :: BamRec -> Maybe [MdOp]
getMd r = case M.lookup (mk_ext_key 'M' 'D') $ b_exts r of
    Just (Text mdfield) -> readMd mdfield
    _                   -> Nothing
  where
    readMd s | L.null s           = return []
             | isDigit (L.head s) = do (n,t) <- maybe (fail "parse error") return $ L.readInt s
                                       (MdNum n :) <$> readMd t
             | L.head s == '^'    = let (a,b) = L.break isDigit (L.tail s)
                                    in (MdDel (L.unpack a) :) <$> readMd b
             | otherwise          = (MdRep (L.head s) :) <$> readMd (L.tail s)


