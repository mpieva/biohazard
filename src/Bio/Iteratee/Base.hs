{-# LANGUAGE TypeFamilies,UndecidableInstances,Rank2Types,ExistentialQuantification #-}

-- |Monadic Iteratees:
-- incremental input parsers, processors and transformers

module Bio.Iteratee.Base (
  -- * Types
  Stream (..)
  ,StreamStatus (..)
  -- ** Exception types
  ,module Bio.Iteratee.Exception
  -- ** Iteratees
  ,Iteratee (..)
  -- * Functions
  -- ** Control functions
  ,run
  ,tryRun
  ,ilift
  ,ifold
  -- ** Creating Iteratees
  ,idone
  ,icont
  ,liftI
  ,idoneM
  ,icontM
  -- ** Stream Functions
  ,setEOF
  -- * Classes
  ,NullPoint(..)
  ,Nullable(..)
)
where

import Bio.Iteratee.Exception
import Bio.Prelude

import Control.Monad.Catch as CIO
import Control.Monad.IO.Class
import Control.Monad.Trans.Class

import qualified Control.Exception    as E
import qualified Data.ByteString      as B
import qualified Data.ByteString.Lazy as L

-- | NullPoint class.  Containers that have a null representation,
-- corresponding to Data.Monoid.mempty.
class NullPoint c where emptyP :: c

instance NullPoint     (Endo a) where emptyP = Endo id
instance NullPoint          [a] where emptyP = []
instance NullPoint B.ByteString where emptyP = B.empty
instance NullPoint L.ByteString where emptyP = L.empty

-- | Nullable container class
class NullPoint c => Nullable c where nullC :: c -> Bool

instance Nullable          [a] where nullC [] = True ; nullC _  = False
instance Nullable B.ByteString where nullC    = B.null
instance Nullable L.ByteString where nullC    = L.null

-- |A stream is a (continuing) sequence of elements bundled in Chunks.
-- The first variant indicates termination of the stream.
-- Chunk a gives the currently available part of the stream.
-- The stream is not terminated yet.
-- The case (null Chunk) signifies a stream with no currently available
-- data but which is still continuing. A stream processor should,
-- informally speaking, ``suspend itself'' and wait for more data
-- to arrive.

data Stream c = EOF (Maybe SomeException) | Chunk c
  deriving (Show, Typeable)

instance (Eq c) => Eq (Stream c) where
  (Chunk c1) == (Chunk c2)           = c1 == c2
  (EOF Nothing) == (EOF Nothing)     = True
  (EOF (Just e1)) == (EOF (Just e2)) = typeOf e1 == typeOf e2
  _ == _                             = False

instance Monoid c => Monoid (Stream c) where
  mempty                        = Chunk mempty
  mappend (EOF mErr)         _  = EOF mErr
  mappend         _  (EOF mErr) = EOF mErr
  mappend (Chunk s1) (Chunk s2) = Chunk (s1 `mappend` s2)

-- |Map a function over a stream.
instance Functor Stream where
  fmap f (Chunk xs) = Chunk $ f xs
  fmap _ (EOF mErr) = EOF mErr

-- |Describe the status of a stream of data.
data StreamStatus =
  DataRemaining
  | EofNoError
  | EofError SomeException
  deriving (Show, Typeable)

-- ----------------------------------------------
-- create exception type hierarchy

-- |Produce the 'EOF' error message.  If the stream was terminated because
-- of an error, keep the error message.
setEOF :: Stream c -> SomeException
setEOF (EOF (Just e)) = e
setEOF _              = toException EofException

-- ----------------------------------------------
-- | Monadic iteratee
newtype Iteratee s m a = Iteratee{ runIter :: forall r.
          (a -> Stream s -> m r) ->
          ((Stream s -> Iteratee s m a) -> Maybe SomeException -> m r) ->
          m r}

-- ----------------------------------------------

idone :: a -> Stream s -> Iteratee s m a
idone a s = Iteratee $ \onDone _ -> onDone a s

icont :: (Stream s -> Iteratee s m a) -> Maybe SomeException -> Iteratee s m a
icont k e = Iteratee $ \_ onCont -> onCont k e

liftI :: (Stream s -> Iteratee s m a) -> Iteratee s m a
liftI k = Iteratee $ \_ onCont -> onCont k Nothing

-- Monadic versions, frequently used by enumerators
idoneM :: Monad m => a -> Stream s -> m (Iteratee s m a)
idoneM x str = return $ Iteratee $ \onDone _ -> onDone x str

icontM
  :: Monad m =>
     (Stream s -> Iteratee s m a)
     -> Maybe SomeException
     -> m (Iteratee s m a)
icontM k e = return $ Iteratee $ \_ onCont -> onCont k e

instance (Functor m) => Functor (Iteratee s m) where
  fmap f m = Iteratee $ \onDone onCont ->
    let od = onDone . f
        oc = onCont . (fmap f .)
    in runIter m od oc

instance (Functor m, Monad m, Nullable s) => Applicative (Iteratee s m) where
    pure x  = idone x (Chunk emptyP)
    {-# INLINE (<*>) #-}
    m <*> a = m >>= flip fmap a

instance (Monad m, Nullable s) => Monad (Iteratee s m) where
  {-# INLINE return #-}
  return x = Iteratee $ \onDone _ -> onDone x (Chunk emptyP)
  {-# INLINE (>>=) #-}
  (>>=) = bindIteratee

{-# INLINE bindIteratee #-}
bindIteratee :: Nullable s
    => Iteratee s m a
    -> (a -> Iteratee s m b)
    -> Iteratee s m b
bindIteratee = self
    where
        self m f = Iteratee $ \onDone onCont ->
             let m_done a (Chunk s)
                   | nullC s     = runIter (f a) onDone onCont
                 m_done a stream = runIter (f a) (const . flip onDone stream) f_cont
                   where f_cont k Nothing = runIter (k stream) onDone onCont
                         f_cont k e       = onCont k e
             in runIter m m_done (onCont . (flip self f .))

instance NullPoint s => MonadTrans (Iteratee s) where
  lift m = Iteratee $ \onDone _ -> m >>= flip onDone (Chunk emptyP)

instance (MonadIO m, Nullable s, NullPoint s) => MonadIO (Iteratee s m) where
  liftIO = lift . liftIO

instance (MonadThrow m, Nullable s, NullPoint s) =>
  MonadThrow (Iteratee s m) where
    throwM e    = lift $ CIO.throwM e

instance (MonadCatch m, Nullable s, NullPoint s) =>
  MonadCatch (Iteratee s m) where
    m `catch` f = Iteratee $ \od oc -> runIter m od oc `CIO.catch` (\e -> runIter (f e) od oc)

-- prior to exceptions-0.6, these were part of MonadCatch
instance (MonadMask m, Nullable s, NullPoint s) => MonadMask (Iteratee s m) where
    mask q      = Iteratee $ \od oc -> CIO.mask $ \u -> runIter (q $ ilift u) od oc
    uninterruptibleMask q = Iteratee $ \od oc -> CIO.uninterruptibleMask $ \u -> runIter (q $ ilift u) od oc


-- |Send 'EOF' to the @Iteratee@ and disregard the unconsumed part of the
-- stream.  If the iteratee is in an exception state, that exception is
-- thrown with 'Control.Exception.throw'.  Iteratees that do not terminate
-- on @EOF@ will throw 'EofException'.
run :: Monad m => Iteratee s m a -> m a
run iter = runIter iter onDone onCont
 where
   onDone  x _        = return x
   onCont  k Nothing  = runIter (k (EOF Nothing)) onDone onCont'
   onCont  _ (Just e) = E.throw e
   onCont' _ Nothing  = E.throw EofException
   onCont' _ (Just e) = E.throw e

-- |Run an iteratee, returning either the result or the iteratee exception.
-- Note that only internal iteratee exceptions will be returned; exceptions
-- thrown with @Control.Exception.throw@ or @Control.Monad.CatchIO.throw@ will
-- not be returned.
--
-- See 'Data.Iteratee.Exception.IFException' for details.
tryRun :: (Exception e, Monad m) => Iteratee s m a -> m (Either e a)
tryRun iter = runIter iter onDone onCont
  where
    onDone  x _ = return $ Right x
    onCont  k Nothing  = runIter (k (EOF Nothing)) onDone onCont'
    onCont  _ (Just e) = return $ maybeExc e
    onCont' _ Nothing  = return $ maybeExc (toException EofException)
    onCont' _ (Just e) = return $ maybeExc e
    maybeExc e = maybe (Left (E.throw e)) Left (fromException e)

-- | Lift a computation in the inner monad of an iteratee.
--
-- A simple use would be to lift a logger iteratee to a monad stack.
--
-- > logger :: Iteratee String IO ()
-- > logger = mapChunksM_ putStrLn
-- >
-- > loggerG :: MonadIO m => Iteratee String m ()
-- > loggerG = ilift liftIO logger
--
-- A more complex example would involve lifting an iteratee to work with
-- interleaved streams.  See the example at 'Data.Iteratee.ListLike.merge'.
ilift ::
  (Monad m, Monad n)
  => (forall r. m r -> n r)
  -> Iteratee s m a
  -> Iteratee s n a
ilift f i = Iteratee $  \od oc ->
  let onDone a str  = return $ Left (a,str)
      onCont k mErr = return $ Right (ilift f . k, mErr)
  in f (runIter i onDone onCont) >>= either (uncurry od) (uncurry oc)

-- | Lift a computation in the inner monad of an iteratee, while threading
-- through an accumulator.
ifold :: (Monad m, Monad n) => (forall r. m r -> acc -> n (r, acc))
      -> acc -> Iteratee s m a -> Iteratee s n (a, acc)
ifold f acc i = Iteratee $ \ od oc -> do
  (r, acc') <- flip f acc $
    runIter i (curry $ return . Left) (curry $ return . Right)
  either (uncurry (od . flip (,) acc'))
         (uncurry (oc . (ifold f acc .))) r
