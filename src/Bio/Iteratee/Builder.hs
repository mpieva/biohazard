{-# LANGUAGE RecordWildCards, FlexibleContexts, BangPatterns, OverloadedStrings #-}
-- | Buffer builder to assemble Bgzf blocks.  (This will probably be
-- renamed.)  The plan is to serialize stuff (BAM and BCF) into a
-- buffer, then Bgzf chunks from the buffer and reuse it.  This /should/
-- avoid redundant copying and relieve some pressure from the garbage
-- collector.  And I hope to plug a mysterious memory leak that doesn't
-- show up in the profiler.
--
-- Exported functions with @unsafe@ in the name resulting in a type of
-- 'Push' omit the bounds checking.  To use them safely, an appropriate
-- 'ensureBuffer' has to precede them.
--
-- XXX  This may not be the most clever way to do it.  According to the
-- reasoning behind the binary-serialise-cbor package, it would be more
-- clever to have a representation of the things we can 'Push' that's
-- similar to a list, and then a function (an Iteratee?) that consumes
-- the list of tokens and fills a buffer.  

module Bio.Iteratee.Builder where

import Bio.Iteratee hiding ( NullPoint ) 
import Bio.Iteratee.Bgzf
import Bio.Prelude
import Data.NullPoint ( NullPoint(..) )
import Data.Primitive.Addr
import Data.Primitive.ByteArray
import Foreign.Marshal.Alloc
import Foreign.Marshal.Utils
import Foreign.Ptr
import Foreign.Storable
import GHC.Exts

import qualified Data.ByteString            as B
import qualified Data.ByteString.Unsafe     as B
import qualified Data.ByteString.Builder    as B ( Builder, toLazyByteString )
import qualified Data.ByteString.Lazy       as B ( foldrChunks )

-- | The 'MutableByteArray' is garbage collected, so we don't get leaks.
-- Once it has grown to a practical size (and the initial 128k should be
-- very practical), we don't get fragmentation either.  We also avoid
-- copies for the most part, since no intermediate 'ByteString's, either
-- lazy or strict have to be allocated.
data BB = BB { buffer :: {-# UNPACK #-} !(MutableByteArray RealWorld)
             , len    :: {-# UNPACK #-} !Int
             , mark   :: {-# UNPACK #-} !Int
             , mark2  :: {-# UNPACK #-} !Int }

-- This still seems to have considerable overhead.  Don't know if this
-- can be improved by effectively inlining IO and turning the BB into an
-- unboxed tuple.  XXX
newtype Push = Push (BB -> IO BB)

instance Monoid Push where
    {-# INLINE mempty #-}
    mempty                  = Push return
    {-# INLINE mappend #-}
    Push a `mappend` Push b = Push (a >=> b)

instance NullPoint Push where
    empty = Push return


-- | Creates a buffer with initial capacity of ~128k.
newBuffer :: IO BB
newBuffer = newPinnedByteArray 128000 >>= \ar -> return $ BB ar 0 0 0

-- | Ensures a given free space in the buffer by doubling its capacity
-- if necessary.
{-# INLINE ensureBuffer #-}
ensureBuffer :: Int -> Push
ensureBuffer n = Push $ \b -> do
    let sz = sizeofMutableByteArray (buffer b)
    if len b + n < sz
       then return b
       else expandBuffer b

expandBuffer :: BB -> IO BB
expandBuffer b = do let sz = sizeofMutableByteArray (buffer b)
                    arr1 <- newPinnedByteArray (sz+sz)
                    copyMutableByteArray arr1 0 (buffer b) 0 (len b)
                    return $ b { buffer = arr1 }

{-# INLINE unsafePushByte #-}
unsafePushByte :: Word8 -> Push
unsafePushByte w = Push $ \b -> do
    writeByteArray (buffer b) (len b) w
    return $ b { len = len b + 1 }

{-# INLINE pushByte #-}
pushByte :: Word8 -> Push
pushByte b = ensureBuffer 1 <> unsafePushByte b

{-# INLINE unsafePushWord32 #-}
unsafePushWord32 :: Word32 -> Push
unsafePushWord32 w = unsafePushByte (fromIntegral $ w `shiftR`  0)
                  <> unsafePushByte (fromIntegral $ w `shiftR`  8)
                  <> unsafePushByte (fromIntegral $ w `shiftR` 16)
                  <> unsafePushByte (fromIntegral $ w `shiftR` 24)

{-# INLINE unsafePushWord16 #-}
unsafePushWord16 :: Word16 -> Push
unsafePushWord16 w = unsafePushByte (fromIntegral $ w `shiftR`  0)
                  <> unsafePushByte (fromIntegral $ w `shiftR`  8)

{-# INLINE pushWord32 #-}
pushWord32 :: Word32 -> Push
pushWord32 w = ensureBuffer 4 <> unsafePushWord32 w

{-# INLINE pushWord16 #-}
pushWord16 :: Word16 -> Push
pushWord16 w = ensureBuffer 2 <> unsafePushWord16 w

{-# INLINE unsafePushByteString #-}
unsafePushByteString :: B.ByteString -> Push
unsafePushByteString bs = Push $ \b ->
    B.unsafeUseAsCStringLen bs $ \(p,ln) -> do
    case mutableByteArrayContents (buffer b) of
        Addr adr -> copyBytes (Ptr adr `plusPtr` len b) p ln
    return $ b { len = len b + ln }

{-# INLINE pushByteString #-}
pushByteString :: B.ByteString -> Push
pushByteString bs = ensureBuffer (B.length bs) <> unsafePushByteString bs

{-# INLINE unsafePushFloat #-}
unsafePushFloat :: Float -> Push
unsafePushFloat f =
    unsafePushWord32 $ unsafeDupablePerformIO $
    alloca $ \b -> poke (castPtr b) f >> peek b

{-# INLINE pushFloat #-}
pushFloat :: Float -> Push
pushFloat f = ensureBuffer 4 <> unsafePushFloat f

{-# INLINE pushBuilder #-}
pushBuilder :: B.Builder -> Push
pushBuilder = B.foldrChunks ((<>) . pushByteString) mempty . B.toLazyByteString

-- | Sets a mark.  This can later be filled in with a record length
-- (used to create BAM records).
{-# INLINE unsafeSetMark #-}
unsafeSetMark :: Push
unsafeSetMark = Push $ \b -> return $ b { len = len b + 4, mark = len b }

{-# INLINE setMark #-}
setMark :: Push
setMark = ensureBuffer 4 <> unsafeSetMark

-- | Ends a record by filling the length into the field that was
-- previously marked.  Terrible things will happen if this wasn't
-- preceded by a corresponding 'setMark'.
{-# INLINE endRecord #-}
endRecord :: Push
endRecord = Push $ \b -> do
    let !l = len b - mark b - 4
    writeByteArray (buffer b) (mark b + 0) (fromIntegral $ shiftR l  0 :: Word8)
    writeByteArray (buffer b) (mark b + 1) (fromIntegral $ shiftR l  8 :: Word8)
    writeByteArray (buffer b) (mark b + 2) (fromIntegral $ shiftR l 16 :: Word8)
    writeByteArray (buffer b) (mark b + 3) (fromIntegral $ shiftR l 24 :: Word8)
    return b

-- | Ends the first part of a record.  The length is filled in *before*
-- the mark, which is specifically done to support the *two* length
-- fields in BCF.  It also remembers the current position.  Horrible
-- things happen if this isn't preceeded by *two* succesive invocations
-- of 'setMark'.
{-# INLINE endRecordPart1 #-}
endRecordPart1 :: Push
endRecordPart1 = Push $ \b -> do
    let !l = len b - mark b - 4
    writeByteArray (buffer b) (mark b - 4) (fromIntegral $ shiftR l  0 :: Word8)
    writeByteArray (buffer b) (mark b - 3) (fromIntegral $ shiftR l  8 :: Word8)
    writeByteArray (buffer b) (mark b - 2) (fromIntegral $ shiftR l 16 :: Word8)
    writeByteArray (buffer b) (mark b - 1) (fromIntegral $ shiftR l 24 :: Word8)
    return $ b { mark2 = len b }

-- | Ends the second part of a record.  The length is filled in at the
-- mark, but computed from the sencond mark only.  This is specifically
-- done to support the *two* length fields in BCF.  Horrible things
-- happen if this isn't preceeded by *two* succesive invocations of
-- 'setMark' and one of 'endRecordPart1'.
{-# INLINE endRecordPart2 #-}
endRecordPart2 :: Push
endRecordPart2 = Push $ \b -> do
    let !l = len b - mark2 b
    writeByteArray (buffer b) (mark b + 0) (fromIntegral $ shiftR l  0 :: Word8)
    writeByteArray (buffer b) (mark b + 1) (fromIntegral $ shiftR l  8 :: Word8)
    writeByteArray (buffer b) (mark b + 2) (fromIntegral $ shiftR l 16 :: Word8)
    writeByteArray (buffer b) (mark b + 3) (fromIntegral $ shiftR l 24 :: Word8)
    return b


{-# INLINE encodeBgzfWith #-}
encodeBgzfWith :: MonadIO m => Int -> Enumeratee Push B.ByteString m b
encodeBgzfWith lv o = newBuffer `ioBind` \bb -> eneeCheckIfDone (liftI . step bb) o
  where
    step bb k (EOF  mx) = finalFlush bb k mx
    step bb k (Chunk (Push p)) = p bb `ioBind` \bb' -> tryFlush bb' 0 k

    tryFlush bb off k
        | len bb - off < maxBlockSize
            = copyMutableByteArray (buffer bb) 0 (buffer bb) off (len bb - off)
              `ioBind_` liftI (step (bb { len = len bb - off
                                        , mark = mark bb - off `max` 0 }) k)

        | otherwise
            = (case mutableByteArrayContents (buffer bb) of
                            Addr adr -> compressChunk lv (Ptr adr `plusPtr` off) (fromIntegral maxBlockSize))
              `ioBind` eneeCheckIfDone (tryFlush bb (off+maxBlockSize)) . k . Chunk

    finalFlush bb k mx
        | len bb < maxBlockSize
            = (case mutableByteArrayContents (buffer bb) of
                            Addr adr -> compressChunk lv (Ptr adr) (fromIntegral $ len bb))
              `ioBind` eneeCheckIfDone (finalFlush2 mx) . k . Chunk

        | otherwise
            = error "WTF?!  This wasn't supposed to happen."

    finalFlush2 mx k = idone (k $ Chunk bgzfEofMarker) (EOF mx)



