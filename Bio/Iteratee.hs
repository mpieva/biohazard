-- | Basically a reexport of @Data.Iteratee@ less the names that clash
-- with @Prelude@ plus a handful of utilities.

{-# LANGUAGE PatternGuards, BangPatterns #-}
module Bio.Iteratee (
    groupStreamBy,
    groupStreamOn,
    i'getString,
    i'lookAhead,
    i'take,
    mapStreamM,
    filterStream,
    filterStreamM,
    foldStream,
    foldStreamM,
    mapChunksMP,

    I.mapStream,
    I.takeWhileE,
    I.tryHead,
    I.isFinished,
    I.heads,
    I.breakE,

    ($==),
    ListLike,
    MonadIO,
    MonadCatchIO,
    lift,
    liftIO,
    (>=>), (<=<),

    enumAuxFile,
    enumInputs,
    enumDefaultInputs,
    defaultBufSize,

    Ordering'(..),
    mergeSortStreams,

    Enumerator',
    Enumeratee',
    mergeEnums',

    module Data.Iteratee.Binary,
    module Data.Iteratee.Char,
    module Data.Iteratee.IO,
    module Data.Iteratee.Iteratee,
    module Data.Iteratee.Parallel
                    ) where

import Bio.Base ( findAuxFile )
import Control.Concurrent
import Control.Monad
import Control.Monad.CatchIO
import Control.Monad.IO.Class
import Control.Monad.Trans.Class
import Data.Iteratee.Binary
import Data.Iteratee.Char
import Data.Iteratee.IO hiding ( defaultBufSize )
import Data.Iteratee.Iteratee
import Data.Iteratee.Parallel
import Data.Monoid
import Data.ListLike ( ListLike )
import System.IO ( stdin )
import System.Environment ( getArgs )

import qualified Data.ListLike as LL
import qualified Data.ByteString as S
import qualified Data.Iteratee as I


-- | Grouping on @Iteratee@s.  @groupStreamOn proj inner outer@ executes
-- @inner (proj e)@, where @e@ is the first input element, to obtain an
-- @Iteratee i@, then passes elements @e@ to @i@ as long as @proj e@
-- produces the same result.  If @proj e@ changes or the input ends, the
-- pair of @proj e@ and the result of @run i@ is passed to @outer@.  At
-- end of input, the resulting @outer@ is returned.
groupStreamOn :: (Monad m, LL.ListLike l e, Eq t1, NullPoint l, Nullable l)
              => (e -> t1)
              -> (t1 -> m (Iteratee l m t2))
              -> Enumeratee l [(t1, t2)] m a
groupStreamOn proj inner = eneeCheckIfDone (liftI . step)
  where
    step outer   (EOF   mx) = idone (liftI outer) $ EOF mx
    step outer c@(Chunk as)
        | LL.null as = liftI $ step outer
        | otherwise  = let x = proj (LL.head as) 
                       in lift (inner x) >>= \i -> step' x i outer c

    -- We want to feed a @Chunk@ to the inner @Iteratee@, which might be
    -- finished.  In that case, we would want to abort, but we cannot,
    -- since the outer iteration is still going on.  So instead we
    -- discard data we would have fed to the inner @Iteratee@.  (Use of
    -- @enumPure1Chunk@ is not appropriate, it would accumulate the
    -- data, just to have it discarded by the @run@ that eventually
    -- happens.
    
    step' c it outer (Chunk as)
        | LL.null as = liftI $ step' c it outer
        | (l,r) <- LL.span ((==) c . proj) as, not (LL.null l) =
            let od a    _str = idoneM a $ EOF Nothing
                oc k Nothing = return $ k (Chunk l)
                oc k       m = icontM k m
            in lift (runIter it od oc) >>= \it' -> step' c it' outer (Chunk r)

    step' c it outer str = 
        lift (run it) >>= \b -> eneeCheckIfDone (\k -> step k str) . outer $ Chunk [(c,b)]


-- | Grouping on @Iteratee@s.  @i'groupBy cmp inner outer@ executes
-- @inner@ to obtain an @Iteratee i@, then passes elements @e@ to @i@ as
-- long as @cmp e' e@, where @e'@ is some preceeding element, is true.
-- Else, the result of @run i@ is passed to @outer@ and @i'groupBy@
-- restarts.  At end of input, the resulting @outer@ is returned.
groupStreamBy :: (Monad m, LL.ListLike l t, NullPoint l, Nullable l)
              => (t -> t -> Bool)
              -> m (Iteratee l m t2)
              -> Enumeratee l [t2] m a
groupStreamBy cmp inner = eneeCheckIfDone (liftI . step)
  where
    step outer    (EOF   mx) = idone (liftI outer) $ EOF mx
    step outer  c@(Chunk as)
        | LL.null as = liftI $ step outer
        | otherwise  = lift inner >>= \i -> step' (LL.head as) i outer c

    step' c it outer (Chunk as)
        | LL.null as = liftI $ step' c it outer
        | (l,r) <- LL.span (cmp c) as, not (LL.null l) =
            let od a    _str = idoneM a $ EOF Nothing
                oc k Nothing = return $ k (Chunk l)
                oc k       m = icontM k m
            in lift (runIter it od oc) >>= \it' -> step' (LL.head l) it' outer (Chunk r)

    step' _ it outer str = 
        lift (run it) >>= \b -> eneeCheckIfDone (\k -> step k str) . outer $ Chunk [b]


i'take :: (Monad m, Nullable s, ListLike s el) => Int -> Enumeratee s s m a
i'take = I.take


-- | Run an Iteratee, collect the input.  When it finishes, return the
-- result along with *all* input.  Effectively allows lookahead.  Be
-- careful, this will eat memory if the @Iteratee@ doesn't return
-- speedily.
i'lookAhead :: Monoid s => Iteratee s m a -> Iteratee s m a
i'lookAhead = go mempty
  where 
    go acc it = Iteratee $ \od oc -> runIter it (\x _ -> od x (Chunk acc)) (oc . step acc)
    
    -- step acc k (Chunk s) | trace (show ("lookahead", msg, S.length acc, S.length s)) False = undefined
    -- step acc k (EOF _) | trace (show ("lookahead", msg, S.length acc, "EOF")) False = undefined

    step acc k c@(Chunk str) = go (acc `mappend` str) (k c)
    step acc k c@(EOF     _) = Iteratee $ \od1 -> runIter (k c) (\x _ -> od1 x (Chunk acc))
                                      

-- | Collects a string of a given length.  Don't use this for long
-- strings, use @Data.Iteratee.ListLike.take@ instead.
i'getString :: Monad m => Int -> Iteratee S.ByteString m S.ByteString
i'getString 0 = idone S.empty (Chunk S.empty)
i'getString n = liftI $ step [] 0
  where
    step acc l c@(EOF _) = icont (step acc l) (Just $ setEOF c)
    step acc l (Chunk c) | l + S.length c >= n = let r = S.concat . reverse $ S.take (n-l) c : acc
                                                 in idone r (Chunk $ S.drop (n-l) c)
                         | otherwise           = liftI $ step (c:acc) (l + S.length c)

infixl 1 $==
-- | Compose an @Enumerator'@ with an @Enumeratee@, giving a new
-- @Enumerator'@.
($==) :: Monad m => Enumerator' hdr input m (Iteratee output m result)
                 -> Enumeratee      input             output m result 
                 -> Enumerator' hdr                   output m result
($==) enum enee iter = run =<< enum (\hdr -> enee $ iter hdr)

-- | Merge two @Enumerator'@s into one.  The header provided by the
-- inner @Enumerator'@ is passed to the output iterator, the header
-- provided by the outer @Enumerator'@ is passed to the merging iteratee
--
-- XXX  Something about those headers is unsatisfactory... there should
--      be an unobtrusive way to combine headers.

mergeEnums' :: (Nullable s2, Nullable s1, Monad m)
            => Enumerator' hi s1 m a                            -- ^ inner enumerator
            -> Enumerator' ho s2 (Iteratee s1 m) a              -- ^ outer enumerator
            -> (ho -> Enumeratee  s2 s1 (Iteratee s1 m) a)      -- ^ merging enumeratee
            -> Enumerator' hi s1 m a
mergeEnums' e1 e2 etee i = e1 $ \hi -> e2 (\ho -> joinI . etee ho $ ilift lift (i hi)) >>= run

-- | Apply a filter predicate to an @Iteratee@.
filterStream :: (Monad m, ListLike s a, NullPoint s) => (a -> Bool) -> Enumeratee s s m r
filterStream = mapChunks . LL.filter

-- | Apply a monadic filter predicate to an @Iteratee@.
filterStreamM :: (Monad m, ListLike s a, Nullable s, NullPoint s) => (a -> m Bool) -> Enumeratee s s m r
filterStreamM k = mapChunksM (go id)
  where
    go acc s | LL.null s = return $! acc LL.empty
             | otherwise = do p <- k (LL.head s)
                              let acc' = if p then LL.cons (LL.head s) . acc else acc
                              go acc' (LL.tail s)

-- | Map a monadic function over an @Iteratee@.
mapStreamM :: (Monad m, ListLike (s el) el, ListLike (s el') el', NullPoint (s el), Nullable (s el), LooseMap s el el')
           => (el -> m el') -> Enumeratee (s el) (s el') m a
mapStreamM = mapChunksM . LL.mapM

-- | Fold a monadic function over an @Iteratee@.
foldStreamM :: (Monad m, Nullable s, ListLike s a) => (b -> a -> m b) -> b -> Iteratee s m b
foldStreamM k = foldChunksM go
  where
    go b s | LL.null s = return b
           | otherwise = k b (LL.head s) >>= \b' -> go b' (LL.tail s)

-- | Fold a function over an @Iteratee@.
foldStream :: (Monad m, Nullable s, ListLike s a) => (b -> a -> b) -> b -> Iteratee s m b
foldStream f = foldChunksM (\b s -> return $! LL.foldl' f b s)

type Enumerator' h eo m b = (h -> Iteratee eo m b) -> m (Iteratee eo m b)
type Enumeratee' h ei eo m b = (h -> Iteratee eo m b) -> Iteratee ei m (Iteratee eo m b)

enumAuxFile :: MonadCatchIO m => FilePath -> Iteratee S.ByteString m a -> m a
enumAuxFile fp it = liftIO (findAuxFile fp) >>= fileDriver it

enumDefaultInputs :: MonadCatchIO m => Enumerator S.ByteString m a
enumDefaultInputs it0 = liftIO getArgs >>= flip enumInputs it0

enumInputs :: MonadCatchIO m => [FilePath] -> Enumerator S.ByteString m a
enumInputs [] = enumHandle defaultBufSize stdin
enumInputs xs = go xs
  where go ("-":fs) = enumHandle defaultBufSize stdin >=> go fs
        go ( f :fs) = enumFile defaultBufSize f >=> go fs
        go [      ] = return

-- | Default buffer size in elements.  This is 1024 in @Iteratee@, which
-- is obviously too small.  Since we want to merge many files, a read
-- should take more time than a seek.  This sets the sensible buffer
-- size to more than about one MB.
defaultBufSize :: Int
defaultBufSize = 2*1024*1024


data Ordering' a = Less | Equal a | NotLess

mergeSortStreams :: (Monad m, ListLike s a, Nullable s) => (a -> a -> Ordering' a) -> Enumeratee s s (Iteratee s m) b
mergeSortStreams comp = eneeCheckIfDone step 
  where
    step out = I.peek >>= \mx -> lift I.peek >>= \my -> case (mx, my) of
        (Just x, Just y) -> case x `comp` y of
            Less    -> do I.drop 1 ;                   eneeCheckIfDone step . out . Chunk $ LL.singleton x
            NotLess -> do            lift (I.drop 1) ; eneeCheckIfDone step . out . Chunk $ LL.singleton y
            Equal z -> do I.drop 1 ; lift (I.drop 1) ; eneeCheckIfDone step . out . Chunk $ LL.singleton z

        (Just  x, Nothing) -> do       I.drop 1  ; eneeCheckIfDone step . out . Chunk $ LL.singleton x
        (Nothing, Just  y) -> do lift (I.drop 1) ; eneeCheckIfDone step . out . Chunk $ LL.singleton y
        (Nothing, Nothing) -> do idone (liftI out) $ EOF Nothing

-- | Map a function over chunks, running in a separate (light weight)
-- thread for each chunk.  MonadIO is needed for the forking.
newtype Ch a = Ch (MVar (Maybe (a, Ch a)))

mapChunksMP :: (MonadIO m, Nullable a) => (a -> IO b) -> Enumeratee a b m c
mapChunksMP f it = do chan <- liftIO $ Ch `liftM` newEmptyMVar
                      eneeCheckIfDone (liftI . go 0 chan chan) it
  where
    maxqueue = 32 :: Int -- arbitrary

    go !_num !chan (Ch !back) k (EOF mx) = do
        -- end the channel, then empty it
        liftIO $ putMVar back Nothing
        let loop (Ch c) k' = do mr <- liftIO $ takeMVar c
                                case mr of 
                                    -- end marker
                                    Nothing -> idone (liftI k') (EOF mx)
                                    Just (r,c') -> eneeCheckIfDone (loop c') . k' $ Chunk r
        loop chan k

    go !num (Ch !chan) (Ch !back) k (Chunk c) = do
        -- First try to get the next finished chunk.  If the queue is
        -- full, don't try, but wait.  If we get something, pass it on
        -- and recurse immediately.
        mnext <- liftIO $ if num >= maxqueue then Just `liftM` takeMVar chan else tryTakeMVar chan
        case mnext of 
            Just (Just (r,chan')) -> eneeCheckIfDone (\k' -> go (num-1) chan' (Ch back) k' (Chunk c)) . k $ Chunk r

            -- end marker... shouldn't actually happen
            Just Nothing -> idone (liftI k) (Chunk c)

            -- if we got nothing, fork and recurse
            Nothing -> do back' <- Ch `liftM` liftIO newEmptyMVar
                          _ <- liftIO . forkIO $ do r <- f c ; putMVar back (Just (r, back'))
                          liftI $ go (num+1) (Ch chan) back' k

