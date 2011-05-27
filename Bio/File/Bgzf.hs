{-# LANGUAGE ForeignFunctionInterface #-}
module Bio.File.Bgzf where

import qualified Data.ByteString as B
import Data.ByteString ( ByteString )
import qualified Data.Iteratee.ListLike as I
import Data.Iteratee.Binary
import Data.Bits
import Control.Monad.IO.Class
import Control.Monad
import Foreign.Marshal.Alloc
import Foreign.Marshal.Utils
import Foreign.Storable
import qualified Data.ByteString.Unsafe as B
import Foreign.C.Types
import Foreign.Ptr
import Data.Word
import Control.Monad.Trans.Class
import Data.Monoid
import Data.Int
import Data.Nullable
import Data.Iteratee.IO
import Data.Iteratee.Char
import Codec.Compression.GZip
import qualified Data.ByteString.Lazy as L

test = enumFile 16300 "/mnt/454/Altaiensis/bwa/catalog/EPO/combined_SNC_anno.tsv.bgz" $
       decompressBgzf $ printLinesUnterminated

decompressBgzf :: MonadIO m => I.Iteratee ByteString m a -> I.Iteratee ByteString m a
decompressBgzf inner = collectI get_bgzf_header $ \csize -> do
                       comp <- get_block $ fromIntegral csize +1
                       let c = B.concat . L.toChunks . decompress $ L.fromChunks [comp]
                       liftIO $ print $ B.length c


                       -- c <- liftIO $ do pu <- mallocBytes (fromIntegral usize)
                                        -- B.unsafeUseAsCStringLen comp $ \(pc, lc) -> do
                                            -- guard (lc == fromIntegral csize + 1)
                                            -- with (fromIntegral usize) $ \plu -> do
                                                -- rc <- zlib_uncompress pu plu pc (fromIntegral lc)
                                                -- unless (rc == 0) . fail $ "fuck, Zlib don't like me: " ++ show rc
                                                -- print (plu, usize)
                                                -- peek plu >>= guard . (fromIntegral usize ==)
                                                -- B.packCStringLen (pu,fromIntegral usize)
                       lift (I.enumPure1Chunk c inner) >>= decompressBgzf

-- Run an Iteratee, collect the input.  When it finishes, apply the
-- continuation and run the result on all input.
collectI :: (Monad m, Monoid s, Nullable s) => I.Iteratee s m a -> (a -> I.Iteratee s m b) -> I.Iteratee s m b
collectI it0 k = go mempty it0
  where 
    go acc it = I.Iteratee $ \od oc -> I.runIter it (onDone od oc) (onCont od oc)
      where
        onDone od oc x _ = I.enumPure1Chunk acc (k x) >>= \it2 ->
                           I.runIter it2 od oc
        onCont od oc k mErr = I.runIter (I.icont (step k) mErr) od oc
        step k c@(I.EOF _) = go acc (k c)
        step k c@(I.Chunk str)
          | I.nullC str            = I.liftI (step k)
          | otherwise            = go (acc `mappend` str) (k c)

get_bgzf_header :: Monad m => I.Iteratee ByteString m Word16
get_bgzf_header = do 31 <- I.head
                     139 <- I.head
                     cm <- I.head
                     flg <- I.head
                     junk <- if flg `testBit` 2 then I.drop 6 else fail "no BGZF"
                     xlen <- endianRead2 LSB 
                     I.joinI $ I.take (fromIntegral xlen) get_bsize

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

get_bsize :: Monad m => I.Iteratee ByteString m Word16
get_bsize = do i1 <- I.head
               i2 <- I.head
               len <- endianRead2 LSB
               if i1 == 66 && i2 == 67 && len == 2 
                  then endianRead2 LSB
                  else I.drop (fromIntegral len) >> get_bsize

get_block :: Monad m => Int -> I.Iteratee ByteString m ByteString
get_block sz = I.liftI $ \s -> case s of
    I.EOF _ -> I.throwErr $ I.setEOF s 
    I.Chunk c | B.length c < sz -> B.append c `liftM` get_block (sz - B.length c)
              | otherwise       -> I.idone (B.take sz c) (I.Chunk (B.drop sz c))


{-
--
--
--
--
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
-}

-- int uncompress (Bytef *dest, uLongf *destLen, const Bytef *source, uLong sourceLen)
-- foreign import ccall unsafe "zlib.h uncompress" zlib_uncompress
    -- :: Ptr CChar -> Ptr CULong -> Ptr CChar -> CULong -> IO CInt
