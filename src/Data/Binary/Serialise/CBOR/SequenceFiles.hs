{-# LANGUAGE ScopedTypeVariables #-}

-- | Support for "CBOR sequence files".
-- Stolen from Duncan Coutts
-- (https://gist.github.com/dcoutts/798812e040a61ad969c27a45549943c0).
--
-- With this, we can write sequences of things into a file, the
-- sequences don't even need to have the same type.  We can read it back
-- as a sequence (now it needs to be homogenous); there should be ways
-- to get single items, too.

module Data.Binary.Serialise.CBOR.SequenceFiles (
   writeBinaryFileSequence,
   appendBinaryFileSequence,
   hPutBinaryFileSequence,
   readBinaryFileSequenceLazy,
   withBinaryFileSequenceLazy,
   withBinaryFileSequence,
 ) where

import           Control.Exception
import           Data.Monoid
import qualified Data.Binary.Serialise.CBOR    as B
import qualified Data.Binary.Serialise.CBOR.Write as B
import qualified Data.ByteString               as BS
import qualified Data.ByteString.Lazy          as LBS
import qualified Data.ByteString.Builder.Extra as BS
import           Data.IORef
import           Foreign
import           GHC.ForeignPtr (mallocPlainForeignPtrBytes)
import           Prelude
import           System.IO
import           System.IO.Error


-- | Write a file consisting of a sequence of items. The items are written out
-- incrementally.
--
-- The file format is simply the binary-encoded items back to back.
--
writeBinaryFileSequence :: B.Serialise a => FilePath -> [a] -> IO ()
writeBinaryFileSequence file xs =
    withFile file WriteMode $ \hnd ->
      hPutBinaryFileSequence hnd xs

-- | Write a file consisting of a sequence of items. The items are written out
-- incrementally.
--
-- The file format is simply the binary-encoded items back to back.
--
appendBinaryFileSequence :: B.Serialise a => FilePath -> [a] -> IO ()
appendBinaryFileSequence file xs =
    withFile file AppendMode $ \hnd ->
      hPutBinaryFileSequence hnd xs


hPutBinaryFileSequence ::  B.Serialise a => Handle -> [a] -> IO ()
hPutBinaryFileSequence hnd xs = do
    buf <- newBuffer BS.defaultChunkSize
    hPutBinaryFileSequenceWithBuffer hnd buf xs


hPutBinaryFileSequenceWithBuffer :: B.Serialise a
                                 => Handle -> Buffer -> [a] -> IO ()
hPutBinaryFileSequenceWithBuffer hnd buf0 =
    go buf0 . BS.runBuilder . B.toBuilder
            . foldr (\x r -> B.encode x <> r) mempty
  where
    go :: Buffer -> BS.BufferWriter -> IO ()
    go buf write = do
      next <- withBuffer buf $ \ptr sz -> do
        -- run the builder, writing into our buffer
        (n, next) <- write ptr sz
        -- so now our buffer contains 'n' bytes
        -- write it all out to the handle leaving our buffer empty
        hPutBuf hnd ptr n
        return next
      case next of
        BS.Done -> return ()
        BS.More minSize write' | bufferSize buf < minSize -> do
          -- very unlikely given our strategy of flushing our buffer every time
          buf' <- newBuffer minSize
          go buf' write'
        BS.More _minSize write' ->
          go buf write'
        BS.Chunk chunk   write' -> do
          BS.hPut hnd chunk
          go buf write'


data Buffer = Buffer {-# UNPACK #-} !(ForeignPtr Word8) {-# UNPACK #-} !Int

bufferSize :: Buffer -> Int
bufferSize (Buffer _fptr len) = len

newBuffer :: Int -> IO Buffer
newBuffer len = do
    fptr <- mallocPlainForeignPtrBytes len
    return $! Buffer fptr len

withBuffer :: Buffer -> (Ptr Word8 -> Int -> IO a) -> IO a
withBuffer (Buffer fptr len) action =
    withForeignPtr fptr $ \ptr -> action ptr len


-- | Read a file consisting of a sequence of items. The items are read in
-- incrementally. The body action is given an action that can be called
-- repeatedly to get each item. It eventually returns @Nothing@.
--
-- The file format is that used by 'writeBinaryFileSequence'.
--
withBinaryFileSequence :: forall a b. B.Serialise a
                       => FilePath
                       -> (IO (Maybe a) -> IO b)
                       -> IO b
withBinaryFileSequence file action = do
    trailingRef <- newIORef BS.empty
    withFile file ReadMode $ \hnd ->
      action (readNextChunk hnd trailingRef)
  where
    readNextChunk :: Handle -> IORef BS.ByteString -> IO (Maybe a)
    readNextChunk hnd trailingRef = do

        initial <- do
          trailing <- readIORef trailingRef
          if BS.null trailing
            then BS.hGetSome hnd BS.defaultChunkSize
            else return trailing

        if BS.null initial
          then return Nothing
          else Just <$> go hnd trailingRef (pushChunk initialDecoder initial)
      where
        initialDecoder :: B.IDecode a
        initialDecoder = B.deserialiseIncremental

    go :: Handle -> IORef BS.ByteString -> B.IDecode a -> IO a
    go hnd trailingRef (B.Partial k) = do
      chunk <- BS.hGetSome hnd BS.defaultChunkSize
      go hnd trailingRef (k (if BS.null chunk then Nothing else Just chunk))

    go _ trailingRef (B.Done trailing _ x) = do
      writeIORef trailingRef trailing
      return x

    go hnd _ (B.Fail _ _ msg) =
      ioError $ mkIOError userErrorType (show msg) (Just hnd) Nothing

-- | Read a file consisting of a sequence of items. The items are read in
-- incrementally, lazily. The body action is given the lazy list of items.
-- These must be consumed within the body action.
--
-- The file format is that used by 'writeBinaryFileSequence'.
--
withBinaryFileSequenceLazy :: forall a b. B.Serialise a
                           => FilePath -> ([a] -> IO b) -> IO b
withBinaryFileSequenceLazy file action = do
    withFile file ReadMode $ \hnd ->
      action . decodeSequence file =<< LBS.hGetContents hnd

-- | Read a file consisting of a sequence of items. The items are read in
-- incrementally, lazily.
--
-- The file format is that used by 'writeBinaryFileSequence'.
--
readBinaryFileSequenceLazy :: forall a. B.Serialise a => FilePath -> IO [a]
readBinaryFileSequenceLazy file =
    decodeSequence file <$> LBS.readFile file

decodeSequence :: forall a. B.Serialise a => FilePath -> LBS.ByteString -> [a]
decodeSequence file =
    go initialDecoder . LBS.toChunks
  where
    initialDecoder :: B.IDecode a
    initialDecoder = B.deserialiseIncremental

    go :: B.IDecode a -> [BS.ByteString] -> [a]
    go (B.Partial k) []             =     go (k Nothing)      []
    go (B.Partial k) (chunk:chunks) =     go (k (Just chunk)) chunks
    go (B.Done trailing _ x) []
                 | BS.null trailing = x : []
    go (B.Done trailing _ x) chunks = x : go initialDecoder  (trailing : chunks)
    go (B.Fail _ _ msg)      _      = throw ioerr
      where ioerr = mkIOError userErrorType (show msg) Nothing (Just file)

pushChunk :: B.IDecode a -> BS.ByteString -> B.IDecode a
pushChunk r inp =
  case r of
    B.Done inp0 p a -> B.Done (inp0 `BS.append` inp) p a
    B.Partial k -> k (Just inp)
    B.Fail inp0 p s -> B.Fail (inp0 `BS.append` inp) p s

