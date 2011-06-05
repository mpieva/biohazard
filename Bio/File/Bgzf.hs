{-# LANGUAGE ForeignFunctionInterface, BangPatterns #-}
module Bio.File.Bgzf where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.Class
-- import Foreign.Marshal.Alloc
-- import Foreign.Marshal.Utils
-- import Foreign.Storable
-- import Foreign.C.Types
-- import Foreign.Ptr
import Data.Bits
import Data.Int                             ( Int64 )
import Data.Iteratee.Binary
import Data.Iteratee.Char                   ( printLinesUnterminated )
import Data.Iteratee.IO
import Data.Monoid
import Data.Word                            ( Word16 )

import qualified Codec.Compression.GZip     as Z
import qualified Data.ByteString            as B
import qualified Data.ByteString.Lazy       as L
-- import qualified Data.ByteString.Unsafe     as B
import qualified Data.Iteratee.ListLike     as I

import Control.Exception.Base

-- one BGZF block: offset and contents
data Block = Block {-# UNPACK #-} !Int64 {-# UNPACK #-} !B.ByteString

instance I.NullPoint Block where empty = Block 0 B.empty
instance I.Nullable Block where nullC (Block _ s) = B.null s

decompress' :: MonadIO m => (I.Iteratee Block m a) -> I.Iteratee B.ByteString m a
decompress' = decompressWith Block 0

decompress :: MonadIO m => (I.Iteratee B.ByteString m a) -> I.Iteratee B.ByteString m a
decompress = decompressWith (\_ s -> s) 0

decompressWith :: MonadIO m => (Int64 -> B.ByteString -> s) -> Int64 -> (I.Iteratee s m a) -> I.Iteratee B.ByteString m a
decompressWith block !off inner = I.isFinished >>= go
  where
    go True = lift $ I.run =<< I.enumEof inner
    go False = collectI get_bgzf_header $ \csize -> do
               comp <- get_block $ fromIntegral csize +1
               -- this is ugly and very roundabout, but works for the time being...
               let c = B.concat . L.toChunks . Z.decompress $ L.fromChunks [comp]
               lift (I.enumPure1Chunk (block off c) inner) >>= decompressWith block (off + fromIntegral csize + 1)

   -- Doesn't work.  Maybe because 'uncompress'
   -- gets confused by the headers?
   -- c <- liftIO $ do pu <- mallocBytes (fromIntegral usize)
                    -- B.unsafeUseAsCStringLen comp $ \(pc, lc) -> do
                        -- guard (lc == fromIntegral csize + 1)
                        -- with (fromIntegral usize) $ \plu -> do
                            -- rc <- zlib_uncompress pu plu pc (fromIntegral lc)
                            -- unless (rc == 0) . fail $ "fuck, Zlib don't like me: " ++ show rc
                            -- print (plu, usize)
                            -- peek plu >>= guard . (fromIntegral usize ==)
                            -- B.packCStringLen (pu,fromIntegral usize)


-- Run an Iteratee, collect the input.  When it finishes, apply the
-- continuation and run the result on all input.
collectI :: (MonadIO m, Monoid s, I.Nullable s) => I.Iteratee s m a -> (a -> I.Iteratee s m b) -> I.Iteratee s m b
collectI it0 k0 = go mempty it0
  where 
    go acc it = I.Iteratee $ \od oc -> I.runIter it (onDone od oc) (onCont od oc)
      where
        onDone od oc x _ = do it2 <- I.enumPure1Chunk acc (k0 x)
                              I.runIter it2 od oc

        onCont od oc k mErr = oc (step k) mErr

        step k c@(I.Chunk str) | I.nullC str = I.liftI (step k)
                               | otherwise   = go (acc `mappend` str) (k c)

        step k c@(I.EOF me) = I.Iteratee $ \od1 oc1 -> I.runIter (k c) (onDone1 me od1 oc1) (onCont1 (step k) oc1)
                                      
        onDone1 me od1 oc1 x _ = do it2 <- I.enumPure1Chunk acc (k0 x)
                                    it2' <- I.enumChunk (I.EOF me) it2
                                    I.runIter it2' od1 oc1

        -- XXX Smells fishy.  If my first Iteratee didn't
        -- produce anything after being passed EOF, what
        -- continuation do I return?  And can it ever be called?
        onCont1 step oc1 k Nothing = oc1 step (Just (SomeException I.DivergentException))
        onCont1 step oc1 k e = oc1 step e


get_bgzf_header :: Monad m => I.Iteratee B.ByteString m Word16
get_bgzf_header = do 31 <- I.head
                     139 <- I.head
                     _cm <- I.head
                     flg <- I.head
                     if flg `testBit` 2 then I.drop 6 else fail "no BGZF"
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

get_bsize :: Monad m => I.Iteratee B.ByteString m Word16
get_bsize = do i1 <- I.head
               i2 <- I.head
               len <- endianRead2 LSB
               if i1 == 66 && i2 == 67 && len == 2 
                  then endianRead2 LSB
                  else I.drop (fromIntegral len) >> get_bsize

get_block :: Monad m => Int -> I.Iteratee B.ByteString m B.ByteString
get_block sz = I.liftI $ \s -> case s of
    I.EOF _ -> I.throwErr $ I.setEOF s 
    I.Chunk c | B.length c < sz -> B.append c `liftM` get_block (sz - B.length c)
              | otherwise       -> I.idone (B.take sz c) (I.Chunk (B.drop sz c))


test = enumFile 16300 "/mnt/ngs_data/101203_SOLEXA-GA04_00007_PEDi_MM_QF_SR/Ibis/BWA/s_5_L3280_sequence_mq_hg19_nohap.bam" $
       decompress' print_block

print_block :: I.Iteratee Block IO ()
print_block = I.liftI step
  where
    step (I.Chunk (Block p s)) = do liftIO $ print $ (p, B.length s)
                                    I.liftI step
    step e@(I.EOF mx) = do liftIO $ putStrLn $ "EOF " ++ show mx
                           I.idone () e

test' = enumFile 16300 "/mnt/454/Altaiensis/bwa/catalog/EPO/combined_SNC_anno.tsv.bgz" $
        decompress printLinesUnterminated

