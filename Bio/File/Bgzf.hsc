{-# LANGUAGE ForeignFunctionInterface, BangPatterns, MultiParamTypeClasses, OverloadedStrings #-}
-- :vim:syn=haskell:

-- | Handling of BGZF files.  Right now, we have an Enumeratee each for
-- input and output.  The input iteratee can optionally supply virtual
-- file offsets, so that seeking is possible.
--
-- Note:  The Zlib bindings are awfully inconvenient for this style.

module Bio.File.Bgzf (
    decompress, decompress', decompressWith, decompressPlain,
    Block(..), compress, maxBlockSize, bgzfEofMarker, 
    liftBlock, getOffset, virtualSeek, isBgzf, isGzip,
    Z.CompressionLevel, Z.noCompression, Z.bestCompression
                     ) where

import Bio.Iteratee
-- import Control.Exception
import Control.Monad
import Foreign.Marshal.Alloc
import Foreign.Storable
import Foreign.C.String
import Foreign.C.Types
import Foreign.Ptr
import Data.Bits
import Data.Monoid
import Data.Word                            ( Word32, Word16, Word8 )

import qualified Codec.Compression.GZip     as Z
import qualified Data.ByteString            as S
import qualified Data.ByteString.Lazy       as L
import qualified Data.ByteString.Unsafe     as S
import qualified Data.Iteratee.ListLike     as I
import qualified Data.ListLike              as LL

#include <zlib.h>

-- | One BGZF block: virtual offset and contents.  Could also be a block
-- of an uncompressed file, if we want to support indexing of
-- uncompressed BAM or some silliness like that.
data Block = Block {-# UNPACK #-} !FileOffset {-# UNPACK #-} !S.ByteString

instance NullPoint Block where empty = Block 0 S.empty
instance Nullable Block where nullC (Block _ s) = S.null s

instance Monoid Block where 
    mempty = empty
    mappend (Block x s) (Block _ t) = Block x (s `S.append` t)
    mconcat [] = empty
    mconcat bs@(Block x _:_) = Block x $ S.concat [s|Block _ s <- bs]

-- | Minimum definition, only needed because @drop@ depends on it
-- indirectly.
instance LL.FoldableLL Block Word8 where
    foldl' f e (Block _ s) = S.foldl' f e s
    foldl f e (Block _ s) = S.foldl f e s
    foldr f e (Block _ s) = S.foldr f e s

-- | Minimum defintion so it works, plus support for @drop@, which was
-- all we really needed...
instance LL.ListLike Block Word8 where
    singleton = Block 0 . S.singleton
    head (Block _ s) = S.head s
    tail (Block o s) = Block (o+1) (S.tail s)
    drop n (Block o s) = Block (o + fromIntegral n) (S.drop n s)
    genericLength (Block _ s) = fromIntegral $ S.length s


-- | Decompresses BGZF into @Block@s.  Each block has a starting offset
-- and is otherwise just a @ByteString@.  Seeking to virtual offsets is
-- supporting iff the underlying stream supports seeking.
decompress' :: Monad m => Enumeratee S.ByteString Block m a
decompress' = decompressWith Block 0

-- | Decompresses BGZF.  The blocks become just @ByteString@s, for
-- consumers who don't want to seek.  Seeking to virtual offsets is
-- supporting iff the underlying stream supports seeking.
decompress :: Monad m => Enumeratee S.ByteString S.ByteString m a
decompress = decompressWith (\_ s -> s) 0

-- | "Decompresses" a plain file.  What's actually happening is that the
-- offset in the input stream is tracked and added to the @ByteString@s
-- giving @Block@s.  This results in the same interface as decompressing
-- actual Bgzf.
decompressPlain :: Monad m => Enumeratee S.ByteString Block m a
decompressPlain = eneeCheckIfDone (liftI . step 0)
  where
    step !o it (Chunk s) = eneeCheckIfDone (liftI . step (o + fromIntegral (S.length s))) . it $ Chunk (Block o s)
    step  _ it (EOF  mx) = idone (liftI it) (EOF mx)

-- | Generic decompression where a function determines how to assemble
-- blocks.
decompressWith :: (Monad m, Monoid s, Nullable s, LL.ListLike s e)
               => ( FileOffset -> S.ByteString -> s )
               -> FileOffset -> Enumeratee S.ByteString s m a
decompressWith blk !off inner = I.isFinished >>= go
  where
    go True = return inner
    go False = do !csize <- maybe (fail "no BGZF") return =<< i'lookAhead get_bgzf_header
                  !comp <- get_block $ fromIntegral csize +1
                  -- this is ugly and very roundabout, but works for the time being...
                  let !c = S.concat . L.toChunks . Z.decompress $ L.fromChunks [comp]
                      !off' = off + fromIntegral csize + 1
                  Iteratee $ \od oc -> do
                      it' <- enumPure1Chunk (blk (off `shiftL` 16) c) inner
                      runIter it' (onDone od) (onCont oc od off')

    -- inner Iteratee was done, so we reconstruct the inner Iteratee and
    -- are done, too.
    onDone od a str = od (idone a str) (Chunk S.empty)

    onCont oc od off' k mx = case mx >>= fromException of
        -- Inner Iteratee continues and either everything is fine or we
        -- don't understand the exception, so we just continue with the
        -- reconstructed inner @Iteratee@.  XXX: Should we propagate the
        -- exception?  How?  (Probably by returning icont and discarding
        -- the current chunk.)
        Nothing -> runIter (decompressWith blk off' (icont k mx)) od oc

        -- inner Iteratee continues and we got a @SeekException@.
        -- This means we issue a seek request, then reenter the
        -- decompression loop.  To seek within the decompressed stream,
        -- we pass a new inner Iteratee that first drops a few bytes,
        -- then calls the original Iteratee.
        Just (SeekException o) -> runIter cont od oc
          where cont = do virtualSeek $ o `shiftR` 16
                          decompressWith blk (o .&. complement 0xffff) $
                              I.drop (fromIntegral $ o .&. 0xffff) >> liftI k

    get_block sz = liftI $ \s -> case s of
        EOF _ -> throwErr $ setEOF s 
        Chunk c | S.length c < sz -> S.append c `liftM` get_block (sz - S.length c)
                | otherwise       -> idone (S.take sz c) (Chunk (S.drop sz c))

    

   -- Doesn't work because 'uncompress' gets confused by the headers
   -- c <- liftIO $ do pu <- mallocBytes (fromIntegral usize)
                    -- S.unsafeUseAsCStringLen comp $ \(pc, lc) -> do
                        -- guard (lc == fromIntegral csize + 1)
                        -- with (fromIntegral usize) $ \plu -> do
                            -- rc <- zlib_uncompress pu plu pc (fromIntegral lc)
                            -- unless (rc == 0) . fail $ "fuck, Zlib don't like me: " ++ show rc
                            -- print (plu, usize)
                            -- peek plu >>= guard . (fromIntegral usize ==)
                            -- S.packCStringLen (pu,fromIntegral usize)
{-
#define Z_OK            0
#define Z_STREAM_END    1
#define Z_NEED_DICT     2
#define Z_ERRNO        (-1)
#define Z_STREAM_ERROR (-2)
#define Z_DATA_ERROR   (-3)
#define Z_MEM_ERROR    (-4)
#define Z_BUF_ERROR    (-5)
#define Z_VERSION_ERROR (-6)
-}


-- | Decodes a BGZF block header and returns the block size if
-- successful.
get_bgzf_header :: Monad m => Iteratee S.ByteString m (Maybe Word16)
get_bgzf_header = do n <- I.heads "\31\139"
                     _cm <- I.head
                     flg <- I.head
                     if flg `testBit` 2 then do
                         I.drop 6
                         xlen <- endianRead2 LSB 
                         it <- I.take (fromIntegral xlen) get_bsize >>= lift . tryRun
                         case it of Left e -> throwErr e
                                    Right s | n == 2 -> return $! Just s
                                    _ -> return Nothing
                      else return Nothing
  where
    get_bsize = do i1 <- I.head
                   i2 <- I.head
                   len <- endianRead2 LSB
                   if i1 == 66 && i2 == 67 && len == 2 
                      then endianRead2 LSB
                      else I.drop (fromIntegral len) >> get_bsize

-- | Seeks in a file that supports it (either plain or BGZF).  
virtualSeek :: Monad m => FileOffset -> Iteratee s m ()
virtualSeek o = icont (idone ()) $ Just $ toException $ SeekException o

-- | Tests whether a stream is in BGZF format.  Does not consume any
-- input.
isBgzf :: Monad m => Iteratee S.ByteString m Bool
isBgzf = liftM check $ checkErr $ i'lookAhead $ get_bgzf_header
  where check (Left         _) = False
        check (Right  Nothing) = False
        check (Right (Just _)) = True

-- | Tests whether a stream is in GZip format.  Also returns @True@ on a
-- Bgzf stream, which is technically a special case of GZip.
isGzip :: Monad m => Iteratee S.ByteString m Bool
isGzip = liftM (either (const False) id) $ checkErr $ i'lookAhead $ test
  where
    test = do n <- I.heads "\31\139"
              I.drop 24
              b <- I.isFinished
              return $ not b && n == 2

-- ------------------------------------------------------------------------- Output

-- | Maximum block size for Bgzf: 64k with some room for headers and
-- uncompressible stuff
maxBlockSize :: Int
maxBlockSize = 65450


-- | The EOF marker for BGZF files.
-- This is just an empty string compressed as BGZF.  Appended to BAM
-- files to indicate their end.
bgzfEofMarker :: S.ByteString
bgzfEofMarker = "\x1f\x8b\x8\x4\0\0\0\0\0\xff\x6\0\x42\x43\x2\0\x1b\0\x3\0\0\0\0\0\0\0\0\0"

-- | Compresses a stream of @ByteString@s into a stream of Bgzf blocks.
-- We accumulate an uncompressed block as long as adding a new chunk to
-- it doesn't exceed the max. block size.  If we receive an empty chunk
-- (used as a flush signal), or if we would exceed the block size, we
-- write out a block.  Then we continue writing until we're below block
-- size.  On EOF, we flush and write the end marker.
--
-- XXX Need a way to write an index "on the side".  Additional output
-- streams?  (Implicitly) pair two @Iteratee@s, similar to @I.pair@?
compress :: MonadIO m => Int -> Enumeratee S.ByteString S.ByteString m a
compress lv = eneeCheckIfDone (liftI . step 0 []) ><> mapChunksMP (liftIO . compress1 lv)
  where
    step    _ acc it c@(EOF _) = step1 it
      where
        step1 i | null acc  = step2 i 
                | otherwise = eneeCheckIfDone step2 . i $ Chunk acc
        step2 i = eneeCheckIfDone step3 . i $ Chunk []
        step3 i = idone (liftI i) c

    step alen acc it (Chunk c) 
        | alen + S.length c < maxBlockSize
            = liftI $ step (alen + S.length c) (c:acc) it

        | S.length c < maxBlockSize 
            = eneeCheckIfDone (liftI . step (S.length c) [c]) . it $ Chunk acc     -- XXX index?

        | otherwise = loop c it -- XXX index?

    loop s i | S.null s  = liftI $ step 0 [] i
             | otherwise = eneeCheckIfDone (loop (S.drop maxBlockSize s)) . i $ Chunk [S.take maxBlockSize s]


-- | Compress a collection of strings into a single BGZF block.

-- Okay, performance was lacking... let's do it again, in a more direct
-- style.  We build out block manually.  First check if the compressed
-- data is going to fit---if not, that's a bug.  Then alloc a buffer,
-- fill with adummy header, alloc a ZStream, compress the pieces we were
-- handed one at a time.  Calculate CRC32, finalize header, construct a
-- byte string, return it.
--
-- We could probably get away with @unsafePerformIO@'ing everything in
-- here, but then again, we only do this when we're writing output
-- anyway.  Hence, run in IO.

compress1 :: Int -> [S.ByteString] -> IO S.ByteString
compress1 _lv [] = return bgzfEofMarker
compress1 lv ss0 = do
    let input_length = sum (map S.length ss0)
    when (input_length > maxBlockSize) $ error "Trying to create too big a BGZF block; this is a bug."
    buf <- mallocBytes 65536

    -- steal header from the EOF marker (length is wrong for now)
    S.unsafeUseAsCString bgzfEofMarker $ \eof ->
        forM_ [0,4..16] $ \o -> do x <- peekByteOff eof o
                                   pokeByteOff buf o (x::Word32)

    -- set up ZStream
    stream <- mallocBytes (#{const sizeof(z_stream)})
    #{poke z_stream, msg}       stream nullPtr
    #{poke z_stream, zalloc}    stream nullPtr
    #{poke z_stream, zfree}     stream nullPtr
    #{poke z_stream, opaque}    stream nullPtr
    #{poke z_stream, next_in}   stream nullPtr
    #{poke z_stream, next_out}  stream (buf `plusPtr` 18)
    #{poke z_stream, avail_in}  stream (0 :: CUInt)
    #{poke z_stream, avail_out} stream (65536-18-8 :: CUInt)
 
    z_check "deflateInit2" $ c_deflateInit2 stream (fromIntegral lv) #{const Z_DEFLATED}
                                            (-15) 8 #{const Z_DEFAULT_STRATEGY}

    -- loop over the fragments.  In reverse order!
    let loop (s:ss) = do 
            crc <- loop ss
            S.unsafeUseAsCString s $ \p ->
              case fromIntegral $ S.length s of
                l | l > 0 -> do
                    #{poke z_stream, next_in} stream p
                    #{poke z_stream, avail_in} stream (l :: CUInt)
                    z_check "deflate" $ c_deflate stream #{const Z_NO_FLUSH}
                    c_crc32 crc p l
                _ -> return crc    
        loop [] = c_crc32 0 nullPtr 0
    crc <- loop ss0
        
    z_check "deflate" $ c_deflate stream #{const Z_FINISH}
    z_check "deflateEnd" $ c_deflateEnd stream

    compressed_length <- (+) (18+8) `fmap` #{peek z_stream, total_out} stream
    when (compressed_length > 65536) $ error "produced too big a block" 
    
    -- set length in header
    pokeByteOff buf 16 (fromIntegral $ (compressed_length-1) .&. 0xff :: Word8)
    pokeByteOff buf 17 (fromIntegral $ (compressed_length-1) `shiftR` 8 :: Word8)

    pokeByteOff buf (compressed_length-8) (fromIntegral crc :: Word32)
    pokeByteOff buf (compressed_length-4) (fromIntegral input_length :: Word32)

    S.unsafePackCStringLen (buf,compressed_length)
  where
    z_check msg code = code >>= \c ->
                       when (c /= #{const Z_OK} && c /= #{const Z_STREAM_END}) $
                       error $ msg ++ " failed: " ++ show c


c_deflateInit2 :: Ptr Word8 -> CInt -> CInt -> CInt -> CInt -> CInt -> IO CInt
c_deflateInit2 z a b c d e = withCAString #{const_str ZLIB_VERSION} $ \versionStr ->
    c_deflateInit2_ z a b c d e versionStr (#{const sizeof(z_stream)} :: CInt)

foreign import ccall unsafe "zlib.h deflateInit2_" c_deflateInit2_ ::
    Ptr Word8 -> CInt -> CInt -> CInt -> CInt -> CInt
		      -> Ptr CChar -> CInt -> IO CInt

foreign import ccall safe "zlib.h deflate" c_deflate ::
    Ptr Word8 -> CInt -> IO CInt

foreign import ccall safe "zlib.h deflateEnd" c_deflateEnd ::
    Ptr Word8 -> IO CInt

foreign import ccall safe "zlib.h crc32" c_crc32 ::
    CULong -> Ptr CChar -> CUInt -> IO CULong

-- ------------------------------------------------------------------------------------------------- utils

-- | Get the current virtual offset.  The virtual address in a BGZF
-- stream contains the offset of the current block in the upper 48 bits
-- and the current offset into that block in the lower 16 bits.  This
-- scheme is compatible with the way BAM files are indexed.
getOffset :: Monad m => Iteratee Block m FileOffset
getOffset = liftI step
  where
    step s@(EOF _) = icont step (Just (setEOF s))
    step s@(Chunk (Block o _)) = idone o s

-- | Runs an @Iteratee@ for @ByteString@s when decompressing BGZF.  Adds
-- internal bookkeeping.
liftBlock :: Monad m => Iteratee S.ByteString m a -> Iteratee Block m a
liftBlock = liftI . step 
  where
    step it (EOF ex) = joinI $ lift $ enumChunk (EOF ex) it
                            
    step it (Chunk (Block l s)) = Iteratee $ \od oc ->
            enumPure1Chunk s it >>= \it' -> runIter it' (onDone od) (oc . step . liftI)
      where
        onDone od hdr (Chunk rest) = od hdr (Chunk $ Block (l + fromIntegral (S.length s-S.length rest)) rest)
        onDone od hdr (EOF     ex) = od hdr (EOF ex)


