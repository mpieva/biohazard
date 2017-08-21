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
  ,enumFdCatch
  ,enumFdRandom
)
where

import Bio.Iteratee.Iteratee
import Bio.Prelude hiding ( bracket )
import Control.Monad.Catch
import Control.Monad.IO.Class
import Data.ByteString.Internal (createAndTrim)
import System.IO (SeekMode(..))

-- | Default buffer size in elements.  This was 1024 in "Data.Iteratee",
-- which is obviously too small.  Since we often want to merge many
-- files, a read should take more time than a seek.  This sets the
-- sensible buffer size to somewhat more than one MB.
defaultBufSize :: Int
defaultBufSize = 2*1024*1024

-- ------------------------------------------------------------------------
-- Binary Random IO enumerators

makefdCallback :: MonadIO m => Int -> Fd -> st -> m (Either SomeException ((Bool, st), Bytes))
makefdCallback bufsize fd st = do
  s <- liftIO . createAndTrim bufsize $ \p ->
       fromIntegral <$> fdReadBuf fd (castPtr p) (fromIntegral bufsize)
  return $ Right ((True, st), s)

-- |The enumerator of a POSIX File Descriptor.  This version enumerates
-- over the entire contents of a file, in order, unless stopped by
-- the iteratee.  In particular, seeking is not supported.
enumFd :: MonadIO m => Int -> Fd -> Enumerator Bytes m a
enumFd bufsize fd = enumFromCallback (makefdCallback bufsize fd) ()

-- |A variant of enumFd that catches exceptions raised by the @Iteratee@.
enumFdCatch
 :: (IException e, MonadIO m)
    => Int
    -> Fd
    -> (e -> m (Maybe EnumException))
    -> Enumerator Bytes m a
enumFdCatch bufsize fd handler = enumFromCallbackCatch (makefdCallback bufsize fd) handler ()

-- |The enumerator of a POSIX File Descriptor: a variation of @enumFd@ that
-- supports RandomIO (seek requests).
enumFdRandom :: MonadIO m => Int -> Fd -> Enumerator Bytes m a
enumFdRandom bs fd iter = enumFdCatch bs fd handler iter
  where
    handler (SeekException off) =
        Nothing <$ (liftIO . fdSeek fd AbsoluteSeek $ fromIntegral off)

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


