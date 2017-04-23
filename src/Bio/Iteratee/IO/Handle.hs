{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |Random and Binary IO with generic Iteratees.  These functions use Handles
-- for IO operations, and are provided for compatibility.  When available,
-- the File Descriptor based functions are preferred as these wastefully
-- allocate memory rather than running in constant space.

module Bio.Iteratee.IO.Handle(
  -- * File enumerators
  enumHandle
  ,enumHandleCatch
  ,enumHandleRandom
  ,enumFile
  ,enumFileRandom
  -- * Iteratee drivers
  ,fileDriverHandle
  ,fileDriverRandomHandle
)

where

import Bio.Iteratee.Iteratee
import Bio.Prelude
import Control.Monad.Catch as CIO
import Control.Monad.IO.Class
import Data.ByteString (packCStringLen)
import Foreign.Marshal.Alloc
import System.IO

-- ------------------------------------------------------------------------
-- Binary Random IO enumerators

makeHandleCallback ::
  MonadIO m =>
  Ptr Word8
  -> Int
  -> Handle
  -> st
  -> m (Either SomeException ((Bool, st), Bytes))
makeHandleCallback p bsize h st = do
  n' <- liftIO (CIO.try $ hGetBuf h p bsize :: IO (Either SomeException Int))
  case n' of
    Left e -> return $ Left e
    Right 0 -> return $ Right ((False, st), emptyP)
    Right n -> liftM (\s -> Right ((True, st), s)) $
                 readFromPtr p (fromIntegral n)
  where
    readFromPtr buf l = liftIO $ packCStringLen (castPtr buf, l)


-- |The (monadic) enumerator of a file Handle.  This version enumerates
-- over the entire contents of a file, in order, unless stopped by
-- the iteratee.  In particular, seeking is not supported.
-- Data is read into a buffer of the specified size.
enumHandle ::
  (MonadIO m, MonadMask m) =>
  Int -- ^Buffer size (number of elements per read)
  -> Handle
  -> Enumerator Bytes m a
enumHandle bufsize h i =
  CIO.bracket (liftIO $ mallocBytes bufsize)
              (liftIO . free)
              (\p -> enumFromCallback (makeHandleCallback p bufsize h) () i)

-- |An enumerator of a file handle that catches exceptions raised by
-- the Iteratee.
enumHandleCatch
 ::
 forall e m a.(IException e,
               MonadIO m, MonadMask m) =>
  Int -- ^Buffer size (number of elements per read)
  -> Handle
  -> (e -> m (Maybe EnumException))
  -> Enumerator Bytes m a
enumHandleCatch bufsize h handler i =
  CIO.bracket (liftIO $ mallocBytes bufsize)
              (liftIO . free)
              (\p -> enumFromCallbackCatch (makeHandleCallback p bufsize h) handler () i)


-- |The enumerator of a Handle: a variation of enumHandle that
-- supports RandomIO (seek requests).
-- Data is read into a buffer of the specified size.
enumHandleRandom ::
 forall m a.(MonadIO m, MonadMask m) =>
  Int -- ^ Buffer size (number of elements per read)
  -> Handle
  -> Enumerator Bytes m a
enumHandleRandom bs h i = enumHandleCatch bs h handler i
  where
    handler (SeekException off) =
       liftM (either
              (Just . EnumException :: IOException -> Maybe EnumException)
              (const Nothing))
             . liftIO . CIO.try $ hSeek h AbsoluteSeek $ fromIntegral off

-- ----------------------------------------------
-- File Driver wrapper functions.

enumFile' :: (MonadIO m, MonadMask m) =>
  (Int -> Handle -> Enumerator s m a)
  -> Int -- ^Buffer size
  -> FilePath
  -> Enumerator s m a
enumFile' enumf bufsize filepath iter = CIO.bracket
  (liftIO $ openBinaryFile filepath ReadMode)
  (liftIO . hClose)
  (flip (enumf bufsize) iter)

enumFile ::
  (MonadIO m, MonadMask m)
  => Int                 -- ^Buffer size
  -> FilePath
  -> Enumerator Bytes m a
enumFile = enumFile' enumHandle

enumFileRandom ::
  (MonadIO m, MonadMask m)
  => Int                 -- ^Buffer size
  -> FilePath
  -> Enumerator Bytes m a
enumFileRandom = enumFile' enumHandleRandom

-- |Process a file using the given @Iteratee@.  This function wraps
-- @enumHandle@ as a convenience.
fileDriverHandle
  :: (MonadIO m, MonadMask m) =>
     Int                      -- ^Buffer size (number of elements)
     -> Iteratee Bytes m a
     -> FilePath
     -> m a
fileDriverHandle bufsize iter filepath =
  enumFile bufsize filepath iter >>= run

-- |A version of @fileDriverHandle@ that supports seeking.
fileDriverRandomHandle
  :: (MonadIO m, MonadMask m) =>
     Int                      -- ^ Buffer size (number of elements)
     -> Iteratee Bytes m a
     -> FilePath
     -> m a
fileDriverRandomHandle bufsize iter filepath =
  enumFileRandom bufsize filepath iter >>= run

