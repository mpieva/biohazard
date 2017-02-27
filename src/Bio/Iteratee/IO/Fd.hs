{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ConstraintKinds #-}

-- |Random and Binary IO with generic Iteratees, using File Descriptors for IO.
-- when available, these are the preferred functions for performing IO as they
-- run in constant space and function properly with sockets, pipes, etc.

module Bio.Iteratee.IO.Fd(
  -- * File enumerators
  -- ** FileDescriptor based enumerators for monadic iteratees
  enumFd
  ,enumFdCatch
  ,enumFdRandom
  ,enumFile
  ,enumFileRandom
  -- * Iteratee drivers
  ,fileDriverFd
  ,fileDriverRandomFd
)

where

import Bio.Iteratee.Binary ()
import Bio.Iteratee.IO.Base
import Bio.Iteratee.Iteratee
import Bio.Iteratee.ReadableChunk
import Bio.Prelude
import Control.Concurrent (yield)
import Control.Monad.Catch as CIO
import Control.Monad.IO.Class
import Foreign.Marshal.Alloc
import Foreign.Ptr
import Foreign.Storable
import System.IO (SeekMode(..))


-- ------------------------------------------------------------------------
-- Binary Random IO enumerators

makefdCallback ::
  (MonadIO m, NullPoint s, ReadableChunk s el) =>
  Ptr el
  -> ByteCount
  -> Fd
  -> st
  -> m (Either SomeException ((Bool, st), s))
makefdCallback p bufsize fd st = do
  n <- liftIO $ myfdRead fd (castPtr p) bufsize
  case n of
    Left  _  -> return $ Left (error "myfdRead failed")
    Right 0  -> liftIO yield >> return (Right ((False, st), emptyP))
    Right n' -> liftM (\s -> Right ((True, st), s)) $
                  readFromPtr p (fromIntegral n')

-- |The enumerator of a POSIX File Descriptor.  This version enumerates
-- over the entire contents of a file, in order, unless stopped by
-- the iteratee.  In particular, seeking is not supported.
enumFd
  :: forall s el m a.(NullPoint s, ReadableChunk s el, MonadIO m, MonadMask m) =>
     Int
     -> Fd
     -> Enumerator s m a
enumFd bs fd iter =
  let bufsize = bs * (sizeOf (undefined :: el))
  in CIO.bracket (liftIO $ mallocBytes bufsize)
                 (liftIO . free)
                 (\p -> enumFromCallback (makefdCallback p (fromIntegral bufsize) fd) () iter)

-- |A variant of enumFd that catches exceptions raised by the @Iteratee@.
enumFdCatch
 :: forall e s el m a.(IException e, NullPoint s, ReadableChunk s el, MonadIO m, MonadMask m)
    => Int
    -> Fd
    -> (e -> m (Maybe EnumException))
    -> Enumerator s m a
enumFdCatch bs fd handler iter =
  let bufsize = bs * (sizeOf (undefined :: el))
  in CIO.bracket (liftIO $ mallocBytes bufsize)
                 (liftIO . free)
                 (\p -> enumFromCallbackCatch (makefdCallback p (fromIntegral bufsize) fd) handler () iter)


-- |The enumerator of a POSIX File Descriptor: a variation of @enumFd@ that
-- supports RandomIO (seek requests).
enumFdRandom
 :: forall s el m a.(NullPoint s, ReadableChunk s el, MonadIO m, MonadMask m) =>
    Int
    -> Fd
    -> Enumerator s m a
enumFdRandom bs fd iter = enumFdCatch bs fd handler iter
  where
    handler (SeekException off) =
      liftM (either
             (const . Just $ enStrExc "Error seeking within file descriptor")
             (const Nothing))
            . liftIO . myfdSeek fd AbsoluteSeek $ fromIntegral off

fileDriver
  :: (MonadIO m, MonadMask m) =>
     (Int -> Fd -> Enumerator s m a)
     -> Int
     -> Iteratee s m a
     -> FilePath
     -> m a
fileDriver enumf bufsize iter filepath = CIO.bracket
  (liftIO $ openFd filepath ReadOnly Nothing defaultFileFlags)
  (liftIO . closeFd)
  (run <=< flip (enumf bufsize) iter)

-- |Process a file using the given @Iteratee@.
fileDriverFd
  :: (NullPoint s, MonadIO m, MonadMask m, ReadableChunk s el) =>
     Int -- ^Buffer size (number of elements)
     -> Iteratee s m a
     -> FilePath
     -> m a
fileDriverFd = fileDriver enumFd

-- |A version of fileDriverFd that supports seeking.
fileDriverRandomFd
  :: (NullPoint s, MonadIO m, MonadMask m, ReadableChunk s el) =>
     Int
     -> Iteratee s m a
     -> FilePath
     -> m a
fileDriverRandomFd = fileDriver enumFdRandom

enumFile' :: (MonadIO m, MonadMask m) =>
  (Int -> Fd -> Enumerator s m a)
  -> Int -- ^Buffer size
  -> FilePath
  -> Enumerator s m a
enumFile' enumf bufsize filepath iter = CIO.bracket
  (liftIO $ openFd filepath ReadOnly Nothing defaultFileFlags)
  (liftIO . closeFd)
  (flip (enumf bufsize) iter)

enumFile ::
  (NullPoint s, MonadIO m, MonadMask m, ReadableChunk s el)
  => Int                 -- ^Buffer size
  -> FilePath
  -> Enumerator s m a
enumFile = enumFile' enumFd

enumFileRandom ::
  (NullPoint s, MonadIO m, MonadMask m, ReadableChunk s el)
  => Int                 -- ^Buffer size
  -> FilePath
  -> Enumerator s m a
enumFileRandom = enumFile' enumFdRandom


