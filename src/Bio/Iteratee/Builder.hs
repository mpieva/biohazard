{-# LANGUAGE UnboxedTuples, RecordWildCards, FlexibleContexts, BangPatterns #-}
-- | Buffer builder to assemble Bgzf blocks.  (This will probably be
-- renamed.)  The plan is to serialize stuff (BAM and BCF) into a
-- buffer, then Bgzf chunks from the buffer and reuse it.  This /should/
-- avoid redundant copying and relieve some pressure from the garbage
-- collector.  And I hope to plug a mysterious memory leak that doesn't
-- show up in the profiler.

module Bio.Iteratee.Builder where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Primitive ( RealWorld )
import Data.Bits
import Data.Map.Strict ( foldrWithKey )
import Data.Monoid
import Data.Primitive.Addr
import Data.Primitive.ByteArray
import GHC.Exts
-- import GHC.Int
-- import GHC.IO     ( IO(..) )
-- import GHC.Prim
import GHC.Word

import qualified Data.Vector.Unboxed as V
import qualified Data.ByteString as B
import qualified Data.ByteString.Unsafe as B

import Bio.Base
import Bio.Bam.Header
import Bio.Bam.Rec
import Bio.Iteratee
import Bio.Iteratee.Bgzf

import Foreign.Marshal.Alloc
import Foreign.Marshal.Utils
import Foreign.Storable
import Foreign.Ptr
import Data.Char ( ord )
import System.IO.Unsafe ( unsafePerformIO )

-- | The 'MutableByteArray' should be garbage collected, so we don't
-- get leaks.  Once it has grown to a practical size (and the initial
-- 128k should be very practical), we don't get fragmentation either.
-- We also avoid copies for the most part, since no intermediate
-- 'ByteString's have to be allocated.
data BB = BB { buffer :: {-# UNPACK #-} !(MutableByteArray RealWorld)
             , len    :: {-# UNPACK #-} !Int
             , mark   :: {-# UNPACK #-} !Int }

newtype Push = Push (BB -> IO BB)

instance Monoid Push where
    mempty  = Push return
    Push a `mappend` Push b = Push (a >=> b)

instance NullPoint Push where
    empty = Push return

instance Nullable Push where
    nullC _ = False -- fake :)

-- | Creates a buffer with initial capacity of ~128k.
newBuffer :: IO BB
newBuffer = newPinnedByteArray 128000 >>= \arr -> return $ BB arr 0 0

-- | Ensures a given free space in the buffer by doubling its capacity
-- if necessary.
ensureBuffer :: Int -> Push
ensureBuffer n = Push $ \b ->
    let sz = sizeofMutableByteArray (buffer b)
    in if len b + n < sz
       then return b
       else do arr1 <- newPinnedByteArray (sz+sz)
               copyMutableByteArray arr1 0 (buffer b) 0 (len b)
               return $ b { buffer = arr1 }

unsafePushByte :: Word8 -> Push
unsafePushByte w = Push $ \b -> do
    writeByteArray (buffer b) (len b) w
    return $ b { len = len b + 1 }

pushByte :: Word8 -> Push
pushByte b = ensureBuffer 1 <> unsafePushByte b

unsafePushWord32 :: Word32 -> Push
unsafePushWord32 w = unsafePushByte (fromIntegral $ w `shiftR`  0)
                  <> unsafePushByte (fromIntegral $ w `shiftR`  8)
                  <> unsafePushByte (fromIntegral $ w `shiftR` 16)
                  <> unsafePushByte (fromIntegral $ w `shiftR` 24)

unsafePushWord16 :: Word16 -> Push
unsafePushWord16 w = unsafePushByte (fromIntegral $ w `shiftR`  0)
                  <> unsafePushByte (fromIntegral $ w `shiftR`  8)

pushWord32 :: Word32 -> Push
pushWord32 w = ensureBuffer 4 <> unsafePushWord32 w

pushWord16 :: Word16 -> Push
pushWord16 w = ensureBuffer 2 <> unsafePushWord16 w

unsafePushByteString :: B.ByteString -> Push
unsafePushByteString bs = Push $ \b ->
    B.unsafeUseAsCStringLen bs $ \(p,ln) -> do
    case mutableByteArrayContents (buffer b) of
        Addr adr -> copyBytes (Ptr adr `plusPtr` len b) p ln
    return $ b { len = len b + ln }

-- | Sets a mark.  This can later be filled in with a record length
-- (used to create BAM records).
unsafeSetMark :: Push
unsafeSetMark = Push $ \b -> return $ b { len = len b + 4, mark = len b }

setMark :: Push
setMark = ensureBuffer 4 <> unsafeSetMark

-- | Ends a record by filling the length into the field that was
-- previously marked.  Terrible things will happen if this wasn't
-- preceded by a corresponding 'setMark'.
endRecord :: Push
endRecord = Push $ \b -> do
    let !l = len b - mark b - 4
    writeByteArray (buffer b) (mark b + 0) (fromIntegral $ shiftR l  0 :: Word8)
    writeByteArray (buffer b) (mark b + 1) (fromIntegral $ shiftR l  8 :: Word8)
    writeByteArray (buffer b) (mark b + 2) (fromIntegral $ shiftR l 16 :: Word8)
    writeByteArray (buffer b) (mark b + 3) (fromIntegral $ shiftR l 24 :: Word8)
    return b


pushBam :: BamRec -> Push
pushBam BamRec{..} = mconcat
    [ ensureBuffer minlength
    , unsafeSetMark
    , unsafePushWord32 (unRefseq b_rname)
    , unsafePushWord32 (fromIntegral b_pos)
    , unsafePushByte (fromIntegral $ B.length b_qname)
    , unsafePushByte (unQ b_mapq)
    , unsafePushWord16 (fromIntegral bin)
    , unsafePushWord16 (fromIntegral $ length $ unCigar b_cigar)
    , unsafePushWord16 (fromIntegral b_flag)
    , unsafePushWord32 (fromIntegral $ V.length b_seq)
    , unsafePushWord32 (unRefseq b_mrnm)
    , unsafePushWord32 (fromIntegral b_mpos)
    , unsafePushWord32 (fromIntegral b_isize)
    , unsafePushByteString b_qname
    , unsafePushByte 0
    , foldMap (unsafePushWord32 . encodeCigar) (unCigar b_cigar)
    , pushSeq b_seq
    , unsafePushByteString b_qual
    , foldrWithKey pushExt mempty b_exts
    , endRecord ]
  where
    bin = distinctBin b_pos (cigarToAlnLen b_cigar)
    minlength = 37 + B.length b_qname + 4 * length (unCigar b_cigar) + B.length b_qual + (V.length b_seq + 1) `shiftR` 1
    encodeCigar (op,l) = fromIntegral $ fromEnum op .|. l `shiftL` 4

    pushSeq :: V.Vector Nucleotides -> Push
    pushSeq v = case v V.!? 0 of
                    Nothing -> mempty
                    Just a  -> case v V.!? 1 of
                        Nothing -> unsafePushByte (unNs a `shiftL` 4)
                        Just b  -> unsafePushByte (unNs a `shiftL` 4 .|. unNs b)
                                   <> pushSeq (V.drop 2 v)

    pushExt :: String -> Ext -> Push -> Push
    pushExt [c,d] e k = case e of
        Text t -> common (4 + B.length t) 'Z' $
                  unsafePushByteString t <> unsafePushByte 0

        Bin  t -> common (4 + B.length t) 'H' $
                  unsafePushByteString t <> unsafePushByte 0

        Char c -> common 4 'A' $ unsafePushByte c

        Float f -> common 7 'f' $ unsafePushWord32 (fromIntegral $ fromFloat f)

        Int i   -> case put_some_int (V.singleton i) of
                        (c,op) -> common 7 c (op i)

        IntArr  ia -> case put_some_int ia of
                        (c,op) -> common (4 * V.length ia) 'B' $ unsafePushByte (fromIntegral $ ord c)
                                  <> unsafePushWord32 (fromIntegral $ V.length ia-1)
                                  <> V.foldr ((<>) . op) mempty ia

        FloatArr fa -> common (4 * V.length fa) 'B' $ unsafePushByte (fromIntegral $ ord 'f')
                       <> unsafePushWord32 (fromIntegral $ V.length fa-1)
                       <> V.foldr ((<>) . unsafePushWord32 . fromFloat) mempty fa
      where
        common l z b = ensureBuffer l <> unsafePushByte (fromIntegral $ ord c)
                   <> unsafePushByte (fromIntegral $ ord d)
                   <> unsafePushByte (fromIntegral $ ord z) <> b <> k

        put_some_int :: V.Vector Int -> (Char, Int -> Push)
        put_some_int is
            | V.all (between        0    0xff) is = ('C', unsafePushByte . fromIntegral)
            | V.all (between   (-0x80)   0x7f) is = ('c', unsafePushByte . fromIntegral)
            | V.all (between        0  0xffff) is = ('S', unsafePushWord16 . fromIntegral)
            | V.all (between (-0x8000) 0x7fff) is = ('s', unsafePushWord16 . fromIntegral)
            | V.all                      (> 0) is = ('I', unsafePushWord32 . fromIntegral)
            | otherwise                           = ('i', unsafePushWord32 . fromIntegral)

        between :: Int -> Int -> Int -> Bool
        between l r x = l <= x && x <= r

        fromFloat :: Float -> Word32
        fromFloat float = unsafePerformIO $ alloca $ \buf ->
                          poke (castPtr buf) float >> peek buf

{-# INLINE encodeBgzfWith #-}
encodeBgzfWith :: MonadIO m => Int -> Enumeratee Push B.ByteString m b
encodeBgzfWith lv o = do bb <- liftIO newBuffer
                         eneeCheckIfDone (liftI . step bb) o
  where
    step bb k (EOF  mx) = finalFlush bb k mx
    step bb k (Chunk (Push p)) = liftIO (p bb) >>= \bb' -> tryFlush bb' 0 k

    tryFlush bb off k
        | len bb - off < maxBlockSize
            = do liftIO $ copyMutableByteArray (buffer bb) 0 (buffer bb) off (len bb - off)
                 liftI $ step (bb { len = len bb - off }) k

        | otherwise
            = do out <- liftIO $ case mutableByteArrayContents (buffer bb) of
                            Addr adr -> compressChunk lv (Ptr adr `plusPtr` off) (fromIntegral maxBlockSize)
                 eneeCheckIfDone (tryFlush bb (off+maxBlockSize)) $ k $ Chunk out

    finalFlush bb k mx
        | len bb < maxBlockSize
            = do out <- liftIO $ case mutableByteArrayContents (buffer bb) of
                            Addr adr -> compressChunk lv (Ptr adr) (fromIntegral $ len bb)
                 eneeCheckIfDone (finalFlush2 mx) $ k $ Chunk out

        | otherwise
            = error "WTF?!  This wasn't supposed to happen."

    finalFlush2 mx k = idone (k $ Chunk bgzfEofMarker) (EOF mx)


