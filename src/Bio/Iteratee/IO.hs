-- |Random and Binary IO with generic Iteratees, using File Descriptors for IO.
-- when available, these are the preferred functions for performing IO as they
-- run in constant space and function properly with sockets, pipes, etc.

module Bio.Iteratee.IO(
  -- * Data
  defaultBufSize
  -- * File enumerators
  ,enumFile
  ,enumFileRandom
  -- * FileDescriptor based enumerators for monadic iteratees
  ,enumFd
  ,enumFdRandom
)
where

import Bio.Iteratee.Iteratee
import Bio.Prelude hiding ( bracket, loop )
import Control.Monad.Catch
import Control.Monad.IO.Class
import System.IO (SeekMode(..))

import qualified Data.ByteString as B

-- | Default buffer size in elements.  This was 1024 in "Data.Iteratee",
-- which is obviously too small.  Since we often want to merge many
-- files, a read should take more time than a seek.  This sets the
-- sensible buffer size to somewhat more than one MB.
defaultBufSize :: Int
defaultBufSize = 2*1024*1024

-- |The enumerator of a POSIX File Descriptor.  This version enumerates
-- over the entire contents of a file, in order, unless stopped by
-- the iteratee.  In particular, seeking is not supported.
enumFd :: MonadIO m => Int -> Fd -> Enumerator Bytes m a
enumFd bufsize fd = loop
  where
    loop iter = runIter iter idoneM onCont

    onCont k j@(Just _) = return (icont k j)
    onCont k   Nothing  = do
        s <- liftIO $ fdGet bufsize fd
        if B.null s then return $ liftI k
                    else loop . k $ Chunk s


-- |The enumerator of a POSIX File Descriptor: a variation of @enumFd@ that
-- supports RandomIO (seek requests).
enumFdRandom :: MonadIO m => Int -> Fd -> Enumerator Bytes m a
enumFdRandom bs fd = loop
  where
    loop iter = runIter iter idoneM onCont

    onCont k Nothing  = do
        s <- liftIO $ fdGet bs fd
        if B.null s then return $ liftI k
                    else loop . k $ Chunk s

    onCont k j@(Just e) = case fromException e of
      Just (SeekException off) -> do
                   liftIO . void $ fdSeek fd AbsoluteSeek (fromIntegral off)
                   loop $ liftI k
      Nothing -> return (icont k j)



enumFile' :: (MonadIO m, MonadMask m) =>
  (Int -> Fd -> Enumerator s m a)
  -> Int -- ^Buffer size
  -> FilePath
  -> Enumerator s m a
enumFile' enumf bufsize filepath iter = bracket
  (liftIO $ openFd filepath ReadOnly Nothing defaultFileFlags)
  (liftIO . closeFd)
  (flip (enumf bufsize) iter)

enumFile ::
  (MonadIO m, MonadMask m)
  => Int                 -- ^Buffer size
  -> FilePath
  -> Enumerator Bytes m a
enumFile = enumFile' enumFd

enumFileRandom ::
  (MonadIO m, MonadMask m)
  => Int                 -- ^Buffer size
  -> FilePath
  -> Enumerator Bytes m a
enumFileRandom = enumFile' enumFdRandom


