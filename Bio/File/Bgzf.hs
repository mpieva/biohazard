{-# LANGUAGE ForeignFunctionInterface, BangPatterns, DeriveDataTypeable, MultiParamTypeClasses #-}

-- | Handling of BGZF files.  Right now, we have an Enumeratee each for
-- input and output.  The input iteratee can optionally supply virtual
-- file offsets, so that seeking is possible.
--
-- Note:  The Zlib bindings are awfully inconvenient for this style.

module Bio.File.Bgzf (
    decompress, decompress', decompressWith, Block(..),
    compress, maxBlockSize, bgzfEofMarker, virtualSeek,
    lookAheadI, liftBlock, getOffset
                     ) where

import Bio.Util
import Control.Exception
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
import Data.Iteratee.Base
import Data.Iteratee.Binary
import Data.Iteratee.Char                   ( printLinesUnterminated )
import Data.Iteratee.IO
import Data.Iteratee.Iteratee
import Data.Monoid
import Data.Typeable
import Data.Word                            ( Word16, Word8 )

import qualified Codec.Compression.GZip     as Z
import qualified Data.ByteString            as S
import qualified Data.ByteString.Lazy       as L
-- import qualified Data.ByteString.Unsafe     as S
import qualified Data.Iteratee.ListLike     as I
import qualified Data.ListLike              as LL

-- | One BGZF block: virtual offset and contents.  Could also be a block
-- of an uncompressed file, if we want to support indexing of
-- uncompressed BAM or some silliness like that.
data Block = Block {-# UNPACK #-} !Int64 {-# UNPACK #-} !S.ByteString

instance I.NullPoint Block where empty = Block 0 S.empty
instance I.Nullable Block where nullC (Block _ s) = S.null s

instance Monoid Block where 
    mempty = I.empty
    mappend (Block x s) (Block _ t) = Block x (s `S.append` t)
    mconcat [] = I.empty
    mconcat bs@(Block x _:_) = Block x $ S.concat [s|Block _ s <- bs]

-- Minimum definition, only needed because @drop@ depends on it
-- indirectly.
instance LL.FoldableLL Block Word8 where
    foldl' f e (Block _ s) = S.foldl' f e s
    foldl f e (Block _ s) = S.foldl f e s
    foldr f e (Block _ s) = S.foldr f e s

-- Minimum defintion so it works, plus support for @drop@, which was all
-- we really needed...
instance LL.ListLike Block Word8 where
    singleton = Block 0 . S.singleton
    head (Block _ s) = S.head s
    tail (Block o s) = Block (o+1) (S.tail s)
    drop n (Block o s) = Block (o + fromIntegral n) (S.drop n s)
    genericLength (Block _ s) = fromIntegral $ S.length s


-- | Decompresses BGZF into @Block@s.  Each block has a starting offset
-- and is otherwise just a @ByteString@.
decompress' :: Monad m => Enumeratee S.ByteString Block m a
decompress' = decompressWith Block 0

-- | Decompresses BGZF.  The blocks become just @ByteString@s, for
-- consumers who don't want to seek.
decompress :: Monad m => Enumeratee S.ByteString S.ByteString m a
decompress = decompressWith (\_ s -> s) 0

decompressWith :: (Monad m, Monoid s, Nullable s, LL.ListLike s e) => (Int64 -> S.ByteString -> s) -> Int64 -> Enumeratee S.ByteString s m a
decompressWith blk !off inner = I.isFinished >>= go
  where
    go True = return inner
    go False = do !csize <- lookAheadI get_bgzf_header
                  !comp <- get_block $ fromIntegral csize +1
                  -- this is ugly and very roundabout, but works for the time being...
                  let !c = S.concat . L.toChunks . Z.decompress $ L.fromChunks [comp]
                      !off' = off + fromIntegral csize + 1
                  Iteratee $ \od oc -> do
                      it' <- enumPure1Chunk (blk (off `shiftL` 16) c) inner
                      runIter it' (onDone od) (onCont oc od off')

    -- inner Iteratee is done, so we reconstruct the inner Iteratee and
    -- are done, too.
    onDone od a str = od (idone a str) (Chunk S.empty)

    onCont oc od off' k mx = case mx >>= fromException of
        -- Inner Iteratee continues and either everything is fine or we
        -- don't understand the exception, so we just continue with the
        -- reconstructed inner @Iteratee@.  XXX: Should we propagate the
        -- exception?  How?
        Nothing -> runIter (decompressWith blk off' (icont k mx)) od oc

        -- inner Iteratee continues and we got a @VirtSeekException@.
        -- This means we issue a seek request, then reenter the
        -- decompression loop.  To seek within the decompressed stream,
        -- we pass a new inner Iteratee that first drops a few bytes,
        -- then calls the original Iteratee.
        Just (VirtSeekException o) -> runIter cont od oc
          where cont = do seek . fromIntegral $ o `shiftR` 16
                          decompressWith blk (o .&. complement 0xffff) $
                              I.drop (fromIntegral $ o .&. 0xffff) >> liftI k
    

   -- Doesn't work.  Maybe because 'uncompress'
   -- gets confused by the headers?
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


get_bgzf_header :: Monad m => Iteratee S.ByteString m Word16
get_bgzf_header = do 31 <- I.head
                     139 <- I.head
                     _cm <- I.head
                     flg <- I.head
                     if flg `testBit` 2 then I.drop 6 else fail "no BGZF"
                     xlen <- endianRead2 LSB 
                     joinI $ I.take (fromIntegral xlen) get_bsize

get_bsize :: Monad m => Iteratee S.ByteString m Word16
get_bsize = do i1 <- I.head
               i2 <- I.head
               len <- endianRead2 LSB
               if i1 == 66 && i2 == 67 && len == 2 
                  then endianRead2 LSB
                  else I.drop (fromIntegral len) >> get_bsize

get_block :: Monad m => Int -> Iteratee S.ByteString m S.ByteString
get_block sz = liftI $ \s -> case s of
    EOF _ -> throwErr $ setEOF s 
    Chunk c | S.length c < sz -> S.append c `liftM` get_block (sz - S.length c)
            | otherwise       -> idone (S.take sz c) (Chunk (S.drop sz c))


-- ------------------------------------------------------------------------- Output

-- | Maximum block size for Bgzf: 64k with some room for headers and
-- uncompressible stuff
maxBlockSize :: Int
maxBlockSize = 65450


-- | The EOF marker for BGZF files.
-- This is just an empty string compressed as BGZF.  Appended to BAM
-- files to indicate their end.
bgzfEofMarker :: S.ByteString
bgzfEofMarker = compress1 []


-- | Compress a stream of @ByteString@s into a stream of Bgzf blocks.
-- We accumulate an uncompressed block as long as adding a new chunk to
-- it doesn't exceed the max. block size.  If we receive an empty chunk
-- (used as a flush signal), or if we would exceed the block size, we
-- write out a block.  Then we continue writing until we're below block
-- size.  On EOF, we flush and write the end marker.
--
-- XXX Need a way to write an index "on the side".  Additional output
-- streams?  (Implicitly) pair two @Iteratee@s, similar to @I.pair@?
compress :: Monad m => Enumeratee S.ByteString S.ByteString m a
compress it0 = icont (step it0 0 []) Nothing
  where
    step it    _ acc c@(EOF _) = lift $ enumPure1Chunk (compress1 acc) it >>= \it' ->       -- XXX index?
                                        enumPure1Chunk bgzfEofMarker it' >>= \it'' ->
                                        enumChunk c it''
    step it alen acc (Chunk c) 
        | alen + S.length c < maxBlockSize
            = icont (step it (alen + S.length c) (c:acc)) Nothing

        | S.length c < maxBlockSize
            = do it' <- lift $ enumPure1Chunk (compress1 acc) it -- XXX index?
                 icont (step it' (S.length c) [c]) Nothing

        | otherwise 
            = do let loop i s | S.null s = return i
                     loop i s = do i' <- lift $ enumPure1Chunk (compress1 [S.take maxBlockSize s]) i
                                   loop i' (S.drop maxBlockSize s)
                 -- XXX index?
                 it' <- loop it c
                 icont (step it' 0 []) Nothing


-- | Compress a single string into a BGZF block.
compress1 :: [S.ByteString] -> S.ByteString
compress1 ss | sum (map S.length ss) > maxBlockSize = error "Don't do that!"
compress1 ss = S.concat (L.toChunks hdr) `S.append` rest
  where
    z = S.concat $ L.toChunks $ Z.compress (L.fromChunks (reverse ss))
    (Right hdr, rest) = runGet patch_header z
    patch_header = do !k <- getWord16le
                      !m <- getWord8
                      !f <- getWord8
                      !t <- getWord32le
                      !xf <- getWord8
                      !_os <- getWord8
                      !xlen <- if f `testBit` 2 then getWord16le else return 0

                      return $! runPut $ do 
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
                            putWord16le . fromIntegral $ S.length z + 5 + 
                                if f `testBit` 2 then 0 else 2

newtype VirtSeekException = VirtSeekException Int64 deriving (Typeable, Show)

instance Exception VirtSeekException where
  toException   = iterExceptionToException
  fromException = iterExceptionFromException

virtualSeek :: ( NullPoint s, Monad m ) => Int64 -> Iteratee s m ()
virtualSeek o = throwRecoverableErr (toException $ VirtSeekException o) (idone ())

-- ------------------------------------------------------------------------------------------------- utils

-- | Get the current virtual offset.  The virtual address in a BGZF
-- stream contains the offset of the current block in the upper 48 bits
-- and the current offset into that block in the lower 16 bits.  This
-- scheme is compatible with the way BAM files are indexed.
getOffset :: Monad m => Iteratee Block m Int64
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

-- ------------------------------------------------------------------------------------------------- tests

print_block :: Iteratee Block IO ()
print_block = liftI step
  where
    step (Chunk (Block p s)) = do liftIO $ putStrLn $ "--> " ++ show (p, S.length s)
                                  print_block
    step e@(EOF mx) = do liftIO $ putStrLn $ "EOF " ++ show mx
                         idone () e

test, test' :: IO ()
test = fileDriverRandom (joinI $ decompress' (virtualSeek 0 >> print_block))
       "/mnt/ngs_data/101203_SOLEXA-GA04_00007_PEDi_MM_QF_SR/Ibis/BWA/s_5_L3280_sequence_mq_hg19_nohap.bam" 

test' = fileDriverRandom (joinI $ decompress (virtualSeek 0 >> printLinesUnterminated))
        "/mnt/454/Altaiensis/bwa/catalog/EPO/combined_SNC_anno.tsv.bgz" 

