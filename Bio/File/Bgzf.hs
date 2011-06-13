{-# LANGUAGE ForeignFunctionInterface, BangPatterns #-}
module Bio.File.Bgzf (
    decompress, decompress', decompressWith, Block(..),
    compress, maxBlockSize, bgzfEofMarker,
    lookAheadI
                     ) where

{-
Handling of BGZF files.  For input, we want:
 - an Iteratee transformer interface
 - decompression of BGZF, GZip and no compression
   (in that case: what virtual adresses do we assign?  Do we even try?)
 - support for seeking.  How?

For output, we want:
 - also an Iteratee transformer
 - writing of uncompressed, GZipped, BGZF'ed files
 - maybe automatic creation of some kind of index?

Other things:
 - the Zlib bindings are awfully inconvenient for this style

Right now, there's just the ITeratee transformer to decompress BGZF into
a sequence of blocks and an Iteratee transformer to compress a sequence
of blocks into BGZF format.
-} 

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.Class
-- import Foreign.Marshal.Alloc
-- import Foreign.Marshal.Utils
-- import Foreign.Storable
-- import Foreign.C.Types
-- import Foreign.Ptr
import Data.Binary.Put
import Data.Binary.Strict.Get
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

-- one BGZF block: virtual offset and contents
data Block = Block {-# UNPACK #-} !Int64 {-# UNPACK #-} !B.ByteString

instance I.NullPoint Block where empty = Block 0 B.empty
instance I.Nullable Block where nullC (Block _ s) = B.null s

-- | Decompresses BGZF into @Block@s.  Each block has a starting offset
-- and is otherwise just a @ByteString@.
decompress' :: Monad m => I.Iteratee Block m a -> I.Iteratee B.ByteString m a
decompress' = decompressWith Block 0

-- | Decompresses BGZF.  The blocks become just @ByteString@s, for
-- consumers who don't want to seek.
decompress :: Monad m => I.Iteratee B.ByteString m a -> I.Iteratee B.ByteString m a
decompress = decompressWith (\_ s -> s) 0

decompressWith :: Monad m => (Int64 -> B.ByteString -> s) -> Int64 -> (I.Iteratee s m a) -> I.Iteratee B.ByteString m a
decompressWith block !off inner = I.isFinished >>= go
  where
    go True = lift $ I.run =<< I.enumEof inner
    go False = do csize <- lookAheadI get_bgzf_header
                  comp <- get_block $ fromIntegral csize +1
                  -- this is ugly and very roundabout, but works for the time being...
                  let c = B.concat . L.toChunks . Z.decompress $ L.fromChunks [comp]
                  lift (I.enumPure1Chunk (block (off `shiftL` 16) c) inner)
                        >>= decompressWith block (off + fromIntegral csize + 1)

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


-- | Run an Iteratee, collect the input.  When it finishes, apply the
-- continuation to its result and run the resulting @Iteratee@ on *all*
-- input.  Effectively allows lookahead.
lookAheadI :: Monoid s => I.Iteratee s m a -> I.Iteratee s m a
lookAheadI = go mempty
  where 
    go acc it = I.Iteratee $ \od oc -> I.runIter it (\x _ -> od x (I.Chunk acc)) (oc . step acc)
    
    step acc k c@(I.Chunk str) = go (acc `mappend` str) (k c)
    step acc k c@(I.EOF     _) = I.Iteratee $ \od1 -> I.runIter (k c) (\x _ -> od1 x (I.Chunk acc))
                                      

get_bgzf_header :: Monad m => I.Iteratee B.ByteString m Word16
get_bgzf_header = do 31 <- I.head
                     139 <- I.head
                     _cm <- I.head
                     flg <- I.head
                     if flg `testBit` 2 then I.drop 6 else fail "no BGZF"
                     xlen <- endianRead2 LSB 
                     I.joinI $ I.take (fromIntegral xlen) get_bsize

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


-- ------------------------------------------------------------------------- Output

-- | Maximum block size for Bgzf: 64k with some room for headers and
-- uncompressible stuff
maxBlockSize :: Int64
maxBlockSize = 65450


-- | The EOF marker for BGZF files.
-- This is just an empty string compressed as BGZF.  Appended to BAM
-- files to indicate their end.
bgzfEofMarker :: B.ByteString
bgzfEofMarker = compress1 L.empty


-- | Compress a stream of @ByteString@s into a stream of Bgzf blocks.
-- We accumulate an uncompressed block as long as adding a new chunk to
-- it doesn't exceed the max. block size.  If we receive an empty chunk
-- (used as a flush signal), or if we would exceed the block size, we
-- write out a block.  Then we continue writing until we're below block
-- size.  On EOF, we flush and write the end marker.
--
-- XXX Need a way to write an index "on the side".  Additional output
-- streams?
compress :: Monad m => I.Iteratee B.ByteString m a -> I.Iteratee B.ByteString m a
compress it0 = I.icont (step it0 L.empty) Nothing
  where
    step it acc c@(I.EOF _) = I.joinIM $ I.enumPure1Chunk (compress1 acc) it >>= \it' ->       -- XXX index?
                                         I.enumPure1Chunk bgzfEofMarker it' >>= \it'' ->
                                         I.enumChunk c it''
    step it acc (I.Chunk c) 
        | L.length acc + fromIntegral (B.length c) < maxBlockSize
            = I.icont (step it (acc `L.append` L.fromChunks [c])) Nothing

        | fromIntegral (B.length c) < maxBlockSize
            = do it' <- lift $ I.enumPure1Chunk (compress1 acc) it -- XXX index?
                 I.icont (step it' (L.fromChunks [c])) Nothing

        | otherwise 
            = do let loop i s | L.null s = return i
                     loop i s = do i' <- lift $ I.enumPure1Chunk (compress1 $ L.take maxBlockSize s) i
                                   loop i' (L.drop maxBlockSize s)
                 -- XXX index?
                 it' <- loop it (L.fromChunks [c])
                 I.icont (step it' L.empty) Nothing


-- | Compress a single string into a BGZF block.
compress1 :: L.ByteString -> B.ByteString
compress1 s | L.length s > maxBlockSize = error "Don't do that!"
compress1 l = B.concat (L.toChunks hdr) `B.append` rest
  where
    z = B.concat $ L.toChunks $ Z.compress l
    (Right hdr, rest) = runGet patch_header z
    patch_header = do k <- getWord16le
                      m <- getWord8
                      f <- getWord8
                      t <- getWord32le
                      xf <- getWord8
                      _os <- getWord8
                      xlen <- if f `testBit` 2 then getWord16le else return 0

                      return $ runPut $ do 
                            putWord16le k
                            putWord8 m
                            putWord8 $ f .|. 4
                            putWord32le t
                            putWord8 xf
                            putWord8 0xff   -- unknown OS
                            putWord16le $ xlen + 6
                            putWord8 66
                            putWord8 67
                            putWord16le 2
                            putWord16le . fromIntegral $ B.length z + 5 + 
                                if f `testBit` 2 then 0 else 2

-- ------------------------------------------------------------------------------------------------- tests

print_block :: I.Iteratee Block IO ()
print_block = I.liftI step
  where
    step (I.Chunk (Block p s)) = do liftIO $ putStrLn $ "--> " ++ show (p, B.length s)
                                    print_block
    step e@(I.EOF mx) = do liftIO $ putStrLn $ "EOF " ++ show mx
                           I.idone () e

test, test' :: IO ()
test = fileDriver (decompress' print_block)
       "/mnt/ngs_data/101203_SOLEXA-GA04_00007_PEDi_MM_QF_SR/Ibis/BWA/s_5_L3280_sequence_mq_hg19_nohap.bam" 

test' = fileDriver (decompress printLinesUnterminated)
        "/mnt/454/Altaiensis/bwa/catalog/EPO/combined_SNC_anno.tsv.bgz" 

