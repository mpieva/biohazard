{-# LANGUAGE ConstraintKinds #-}

-- |Random and Binary IO with generic Iteratees.

module Bio.Iteratee.IO(
  -- * Data
  defaultBufSize,
  -- * File enumerators
  -- ** Handle-based enumerators
  H.enumHandle,
  H.enumHandleRandom,
  enumFile,
  enumFileRandom,
  -- ** FileDescriptor based enumerators
  FD.enumFd,
  FD.enumFdRandom,
  -- * Iteratee drivers
  --   These are FileDescriptor-based on POSIX systems, otherwise they are
  --   Handle-based.  The Handle-based drivers are accessible on POSIX systems
  --   at Data.Iteratee.IO.Handle
  fileDriver,
  fileDriverVBuf,
  fileDriverRandom,
  fileDriverRandomVBuf,
)

where

import Bio.Iteratee.ReadableChunk
import Bio.Iteratee.Iteratee
import Bio.Iteratee.Binary ()
import Control.Monad.Catch
import Control.Monad.IO.Class
import Prelude

import qualified Bio.Iteratee.IO.Handle as H
import qualified Bio.Iteratee.IO.Fd as FD


-- | Default buffer size in elements.  This was 1024 in "Data.Iteratee",
-- which is obviously too small.  Since we often want to merge many
-- files, a read should take more time than a seek.  This sets the
-- sensible buffer size to somewhat more than one MB.
defaultBufSize :: Int
defaultBufSize = 2*1024*1024


-- If Posix is available, use the fileDriverRandomFd as fileDriverRandom.  Otherwise, use a handle-based variant.
enumFile
  :: (MonadIO m, MonadMask m, NullPoint s, ReadableChunk s el) =>
     Int
     -> FilePath
     -> Enumerator s m a
enumFile = FD.enumFile

enumFileRandom
  :: (MonadIO m, MonadMask m, NullPoint s, ReadableChunk s el) =>
     Int
     -> FilePath
     -> Enumerator s m a
enumFileRandom = FD.enumFileRandom

-- |Process a file using the given Iteratee.  This function wraps
-- enumFd as a convenience.
fileDriver
  :: (MonadIO m, MonadMask m, NullPoint s, ReadableChunk s el) =>
     Iteratee s m a
     -> FilePath
     -> m a
fileDriver = FD.fileDriverFd defaultBufSize

-- |A version of fileDriver with a user-specified buffer size (in elements).
fileDriverVBuf
  :: (MonadIO m, MonadMask m, NullPoint s, ReadableChunk s el) =>
     Int
     -> Iteratee s m a
     -> FilePath
     -> m a
fileDriverVBuf = FD.fileDriverFd

-- |Process a file using the given Iteratee.  This function wraps
-- enumFdRandom as a convenience.
fileDriverRandom
  :: (MonadIO m, MonadMask m, NullPoint s, ReadableChunk s el) =>
     Iteratee s m a
     -> FilePath
     -> m a
fileDriverRandom = FD.fileDriverRandomFd defaultBufSize

fileDriverRandomVBuf
  :: (MonadIO m, MonadMask m, NullPoint s, ReadableChunk s el) =>
     Int
     -> Iteratee s m a
     -> FilePath
     -> m a
fileDriverRandomVBuf = FD.fileDriverRandomFd

