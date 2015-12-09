{-# LANGUAGE PatternGuards, BangPatterns, DeriveDataTypeable #-}

-- | Basically a reexport of "Data.Iteratee" less the names that clash
-- with "Prelude" plus a handful of utilities.

module Bio.Iteratee (
    groupStreamBy,
    groupStreamOn,
    iGetString,
    iterGet,
    iterLoop,
    iLookAhead,
    headStream,
    peekStream,
    takeStream,
    dropStream,
    mapStreamM,
    mapStreamM_,
    filterStream,
    filterStreamM,
    foldStream,
    foldStreamM,
    zipStreams,
    protectTerm,
    concatMapStream,
    mapMaybeStream,
    parMapChunksIO,
    progressNum,

    I.mapStream,
    I.takeWhileE,
    I.tryHead,
    I.isFinished,
    I.heads,
    I.breakE,

    ($==),
    ioBind, ioBind_,
    ListLike,
    MonadIO, MonadMask,
    lift, liftIO,
    (>=>), (<=<),
    stdin, stdout, stderr,

    enumAuxFile,
    enumInputs,
    enumDefaultInputs,
    defaultBufSize,

    Ordering'(..),
    mergeSortStreams,

    Enumerator',
    Enumeratee',
    mergeEnums',

    QQ(..),
    emptyQ,
    lengthQ,
    pushQ,
    popQ,
    cancelAll,

    ParseError(..),
    parserToIteratee,
    stream2vector,
    stream2vectorN,

    Fd,
    withFileFd,
    module Data.Iteratee.Binary,
    module Data.Iteratee.Char,
    module Data.Iteratee.IO,
    module Data.Iteratee.Iteratee
        ) where

import Bio.Base                             ( findAuxFile )
import Bio.Util                             ( showNum )
import Control.Concurrent.Async             ( Async, async, wait, cancel )
import Control.Monad
import Control.Monad.Catch
import Control.Monad.IO.Class
import Control.Monad.Trans.Class
import Data.Binary.Get
import Data.Iteratee.Binary
import Data.Iteratee.Char
import Data.Iteratee.IO              hiding ( defaultBufSize )
import Data.Iteratee.Iteratee        hiding ( identity )
import Data.ListLike                        ( ListLike )
import Data.Monoid
import Data.Typeable
import System.IO                            ( stdin, stdout, stderr, hIsTerminalDevice )
import System.Environment                   ( getArgs )
import System.Mem                           ( performGC )
import System.Posix                         ( Fd, openFd, closeFd, OpenMode(..), defaultFileFlags )

import qualified Data.Attoparsec.ByteString     as A
import qualified Data.ByteString                as S
import qualified Data.Iteratee                  as I
import qualified Data.ListLike                  as LL
import qualified Data.Vector.Generic            as VG
import qualified Data.Vector.Generic.Mutable    as VM

-- | Grouping on 'Iteratee's.  @groupStreamOn proj inner outer@ executes
-- @inner (proj e)@, where @e@ is the first input element, to obtain an
-- 'Iteratee' @i@, then passes elements @e@ to @i@ as long as @proj e@
-- produces the same result.  If @proj e@ changes or the input ends, the
-- pair of @proj e@ and the result of @run i@ is passed to @outer@.  At
-- end of input, the resulting @outer@ is returned.
groupStreamOn :: (Monad m, LL.ListLike l e, Eq t1, NullPoint l, Nullable l)
              => (e -> t1)
              -> (t1 -> m (Iteratee l m t2))
              -> Enumeratee l [(t1, t2)] m a
groupStreamOn proj inner = eneeCheckIfDonePass (icont . step)
  where
    step outer   (EOF   mx) = idone (liftI outer) $ EOF mx
    step outer c@(Chunk as)
        | LL.null as = liftI $ step outer
        | otherwise  = let x = proj (LL.head as)
                       in lift (inner x) >>= \i -> step' x i outer c

    -- We want to feed a 'Chunk' to the inner 'Iteratee', which might be
    -- finished.  In that case, we would want to abort, but we cannot,
    -- since the outer iteration is still going on.  So instead we
    -- discard data we would have fed to the inner 'Iteratee'.  (Use of
    -- 'enumPure1Chunk' is not appropriate, it would accumulate the
    -- data, just to have it discarded by the 'run' that eventually
    -- happens.

    step' c it outer (Chunk as)
        | LL.null as = liftI $ step' c it outer
        | (l,r) <- LL.span ((==) c . proj) as, not (LL.null l) =
            let od a    _str = idoneM a $ EOF Nothing
                oc k Nothing = return $ k (Chunk l)
                oc k       m = icontM k m
            in lift (runIter it od oc) >>= \it' -> step' c it' outer (Chunk r)

    step' c it outer str =
        lift (run it) >>= \b -> eneeCheckIfDone (`step` str) . outer $ Chunk [(c,b)]


-- | Grouping on 'Iteratee's.  @groupStreamBy cmp inner outer@ executes
-- @inner@ to obtain an 'Iteratee' @i@, then passes elements @e@ to @i@
-- as long as @cmp e0 e@, where @e0@ is some preceeding element, is
-- true.  Else, the result of @run i@ is passed to @outer@ and
-- 'groupStreamBy' restarts.  At end of input, the resulting @outer@ is
-- returned.
groupStreamBy :: (Monad m, LL.ListLike l t, NullPoint l, Nullable l)
              => (t -> t -> Bool)
              -> m (Iteratee l m t2)
              -> Enumeratee l [t2] m a
groupStreamBy cmp inner = eneeCheckIfDonePass (icont . step)
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
        lift (run it) >>= \b -> eneeCheckIfDone (`step` str) . outer $ Chunk [b]


-- | Take a prefix of a stream, the equivalent of 'Data.List.take'.
{-# INLINE takeStream #-}
takeStream :: (Monad m, Nullable s, ListLike s el) => Int -> Enumeratee s s m a
takeStream = I.take

-- | Take first element of a stream or fail.
{-# INLINE headStream #-}
headStream :: ListLike s el => Iteratee s m el
headStream = I.head

{-# INLINE peekStream #-}
peekStream :: ListLike s el => Iteratee s m (Maybe el)
peekStream = I.peek

{-# INLINE dropStream #-}
dropStream :: (Nullable s, ListLike s el) => Int -> Iteratee s m ()
dropStream = I.drop

-- | Run an Iteratee, collect the input.  When it finishes, return the
-- result along with *all* input.  Effectively allows lookahead.  Be
-- careful, this will eat memory if the @Iteratee@ doesn't return
-- speedily.
iLookAhead :: Monoid s => Iteratee s m a -> Iteratee s m a
iLookAhead = go mempty
  where
    go acc it = Iteratee $ \od oc -> runIter it (\x _ -> od x (Chunk acc)) (oc . step acc)

    step acc k c@(Chunk str) = go (acc `mappend` str) (k c)
    step acc k c@(EOF     _) = Iteratee $ \od1 -> runIter (k c) (\x _ -> od1 x (Chunk acc))


-- | Collects a string of a given length.  Don't use this for long
-- strings, use 'takeStream' instead.
iGetString :: Monad m => Int -> Iteratee S.ByteString m S.ByteString
iGetString 0 = idone S.empty (Chunk S.empty)
iGetString n = liftI $ step [] 0
  where
    step acc l c@(EOF _) = icont (step acc l) (Just $ setEOF c)
    step acc l (Chunk c) | l + S.length c >= n = let r = S.concat . reverse $ S.take (n-l) c : acc
                                                 in idone r (Chunk $ S.drop (n-l) c)
                         | otherwise           = liftI $ step (c:acc) (l + S.length c)

-- | Repeatedly apply an 'Iteratee' to a value until end of stream.
-- Returns the final value.
iterLoop :: (Nullable s, Monad m) => (a -> Iteratee s m a) -> a -> Iteratee s m a
iterLoop it a = do e <- I.isFinished
                   if e then return a
                        else it a >>= iterLoop it

-- | Convert a 'Get' into an 'Iteratee'.  The 'Get' is applied once, the
-- decoded data is returned, unneded input remains in the stream.
iterGet :: Monad m => Get a -> Iteratee S.ByteString m a
iterGet = go . runGetIncremental
  where
    go (Fail  _ _ err) = throwErr (iterStrExc err)
    go (Done rest _ a) = idone a (Chunk rest)
    go (Partial   dec) = liftI $ \ck -> case ck of
        Chunk s -> go (dec $ Just s)
        EOF  mx -> case dec Nothing of
            Fail  _ _ err -> throwErr (iterStrExc err)
            Partial     _ -> throwErr (iterStrExc "<partial>")
            Done rest _ a | S.null rest -> idone a (EOF mx)
                          | otherwise   -> idone a (Chunk rest)

{-# INLINE ioBind #-}
-- | Lifts an IO action and combines it with a continution.
-- @ioBind m f@ is the same as @liftIO m >>= f@, but does not create a
-- 'Nullable' constraint on the stream type.
infixl 1 `ioBind`
ioBind :: MonadIO m => IO a -> (a -> Iteratee s m b) -> Iteratee s m b
ioBind m f = Iteratee $ \onDone onCont -> liftIO m >>= \a -> runIter (f a) onDone onCont

{-# INLINE ioBind_ #-}
-- | Lifts an IO action, ignores its result, and combines it with a continution.
-- @ioBind_ m f@ is the same as @liftIO m >> f@, but does not create a
-- 'Nullable' constraint on the stream type.
infixl 1 `ioBind_`
ioBind_ :: MonadIO m => IO a -> Iteratee s m b -> Iteratee s m b
ioBind_ m f = Iteratee $ \onDone onCont -> liftIO m >> runIter f onDone onCont

infixl 1 $==
{-# INLINE ($==) #-}
-- | Compose an 'Enumerator\'' with an 'Enumeratee', giving a new
-- 'Enumerator\''.
($==) :: Monad m => Enumerator' hdr input m (Iteratee output m result)
                 -> Enumeratee      input             output m result
                 -> Enumerator' hdr                   output m result
($==) enum enee iter = run =<< enum (enee . iter)

-- | Merge two 'Enumerator\''s into one.  The header provided by the
-- inner 'Enumerator\'' is passed to the output iterator, the header
-- provided by the outer 'Enumerator\'' is passed to the merging iteratee
--
-- XXX  Something about those headers is unsatisfactory... there should
--      be an unobtrusive way to combine headers.

{-# INLINE mergeEnums' #-}
mergeEnums' :: (Nullable s2, Nullable s1, Monad m)
            => Enumerator' hi s1 m a                            -- ^ inner enumerator
            -> Enumerator' ho s2 (Iteratee s1 m) a              -- ^ outer enumerator
            -> (ho -> Enumeratee  s2 s1 (Iteratee s1 m) a)      -- ^ merging enumeratee
            -> Enumerator' hi s1 m a
mergeEnums' e1 e2 etee i = e1 $ \hi -> e2 (\ho -> joinI . etee ho $ ilift lift (i hi)) >>= run

{-# INLINE concatMapStream #-}
concatMapStream :: (Monad m, ListLike s a, NullPoint s, ListLike t b) => (a -> t) -> Enumeratee s t m r
concatMapStream = mapChunks . LL.concatMap

{-# INLINE mapMaybeStream #-}
mapMaybeStream :: (Monad m, ListLike s a, NullPoint s, ListLike t b) => (a -> Maybe b) -> Enumeratee s t m r
mapMaybeStream f = mapChunks mm
  where
    mm l = if LL.null l then LL.empty else
           case f (LL.head l) of Nothing -> mm (LL.tail l)
                                 Just b  -> LL.cons b $ mm (LL.tail l)

-- | Apply a filter predicate to an 'Iteratee'.
{-# INLINE filterStream #-}
filterStream :: (Monad m, ListLike s a, NullPoint s) => (a -> Bool) -> Enumeratee s s m r
filterStream = mapChunks . LL.filter

-- | Apply a monadic filter predicate to an 'Iteratee'.
{-# INLINE filterStreamM #-}
filterStreamM :: (Monad m, ListLike s a, Nullable s, NullPoint s) => (a -> m Bool) -> Enumeratee s s m r
filterStreamM k = mapChunksM (go id)
  where
    go acc s | LL.null s = return $! acc LL.empty
             | otherwise = do p <- k (LL.head s)
                              let acc' = if p then LL.cons (LL.head s) . acc else acc
                              go acc' (LL.tail s)

-- | Map a monadic function over an 'Iteratee'.
{-# INLINE mapStreamM #-}
mapStreamM :: (Monad m, ListLike (s el) el, ListLike (s el') el', NullPoint (s el), Nullable (s el), LooseMap s el el')
           => (el -> m el') -> Enumeratee (s el) (s el') m a
mapStreamM = mapChunksM . LL.mapM

-- | Map a monadic function over an 'Iteratee', discarding the results.
{-# INLINE mapStreamM_ #-}
mapStreamM_ :: (Monad m, Nullable s, ListLike s el) => (el -> m b) -> Iteratee s m ()
mapStreamM_ = mapChunksM_ . LL.mapM_

-- | Fold a monadic function over an 'Iteratee'.
{-# INLINE foldStreamM #-}
foldStreamM :: (Monad m, Nullable s, ListLike s a) => (b -> a -> m b) -> b -> Iteratee s m b
foldStreamM k = foldChunksM go
  where
    go b s | LL.null s = return b
           | otherwise = k b (LL.head s) >>= \b' -> go b' (LL.tail s)

-- | Fold a function over an 'Iteratee'.
foldStream :: (Monad m, Nullable s, ListLike s a) => (b -> a -> b) -> b -> Iteratee s m b
foldStream f = foldChunksM (\b s -> return $! LL.foldl' f b s)

-- | Apply two 'Iteratee's to the same stream.
zipStreams :: (Nullable s, ListLike s el, Monad m)
           => Iteratee s m a -> Iteratee s m b -> Iteratee s m (a, b)
zipStreams = I.zip

type Enumerator' h eo m b = (h -> Iteratee eo m b) -> m (Iteratee eo m b)
type Enumeratee' h ei eo m b = (h -> Iteratee eo m b) -> Iteratee ei m (Iteratee eo m b)

enumAuxFile :: (MonadIO m, MonadMask m) => FilePath -> Iteratee S.ByteString m a -> m a
enumAuxFile fp it = liftIO (findAuxFile fp) >>= fileDriver it

enumDefaultInputs :: (MonadIO m, MonadMask m) => Enumerator S.ByteString m a
enumDefaultInputs it0 = liftIO getArgs >>= flip enumInputs it0

enumInputs :: (MonadIO m, MonadMask m) => [FilePath] -> Enumerator S.ByteString m a
enumInputs [] = enumHandle defaultBufSize stdin
enumInputs xs = go xs
  where go ("-":fs) = enumHandle defaultBufSize stdin >=> go fs
        go ( f :fs) = enumFile defaultBufSize f >=> go fs
        go [      ] = return

-- | Default buffer size in elements.  This is 1024 in "Data.Iteratee",
-- which is obviously too small.  Since we want to merge many files, a
-- read should take more time than a seek.  This sets the sensible
-- buffer size to more than about one MB.
defaultBufSize :: Int
defaultBufSize = 2*1024*1024


data Ordering' a = Less | Equal a | NotLess

mergeSortStreams :: (Monad m, ListLike s a, Nullable s) => (a -> a -> Ordering' a) -> Enumeratee s s (Iteratee s m) b
mergeSortStreams comp = eneeCheckIfDone step
  where
    step out = peekStream >>= \mx -> lift peekStream >>= \my -> case (mx, my) of
        (Just x, Just y) -> case x `comp` y of
            Less    -> do I.drop 1 ;                   eneeCheckIfDone step . out . Chunk $ LL.singleton x
            NotLess -> do            lift (I.drop 1) ; eneeCheckIfDone step . out . Chunk $ LL.singleton y
            Equal z -> do I.drop 1 ; lift (I.drop 1) ; eneeCheckIfDone step . out . Chunk $ LL.singleton z

        (Just  x, Nothing) -> do       I.drop 1  ; eneeCheckIfDone step . out . Chunk $ LL.singleton x
        (Nothing, Just  y) -> do lift (I.drop 1) ; eneeCheckIfDone step . out . Chunk $ LL.singleton y
        (Nothing, Nothing) -> idone (liftI out) $ EOF Nothing


-- | Parallel map of an IO action over the elements of a stream
--
-- This 'Enumeratee' applies an 'IO' action to every chunk of the input
-- stream.  These 'IO' actions are run asynchronously in a limited
-- parallel way.  Don't forget to `evaluate`

parMapChunksIO :: (MonadIO m, Nullable s) => Int -> (s -> IO t) -> Enumeratee s t m a
parMapChunksIO np f = eneeCheckIfDonePass (go emptyQ)
  where
    -- check if the queue is full
    go !qq k (Just e) = cancelAll qq >> icont (go' emptyQ k) (Just e)
    go !qq k Nothing = case popQ qq of
        Just (a,qq') | lengthQ qq == np -> liftIO (wait a) >>= eneeCheckIfDonePass (go qq') . k . Chunk
        _                               -> liftI $ go' qq k

    -- we have room for input
    go' !qq k (EOF  mx) = do a <- liftIO (async (f empty))
                             goE mx (pushQ a qq) k Nothing
    go' !qq k (Chunk c) = do a <- liftIO (async (f c))
                             go (pushQ a qq) k Nothing

    -- input ended, empty the queue
    goE  _ !qq k (Just e) = cancelAll qq >> icont (go' emptyQ k) (Just e)
    goE mx !qq k Nothing = case popQ qq of
        Nothing      -> idone (liftI k) (EOF mx)
        Just (a,qq') -> liftIO (wait a) >>= eneeCheckIfDonePass (goE mx qq') . k . Chunk

-- | Protects the terminal from binary junk.  If @i@ is an 'Iteratee'
-- that might write binary to 'stdout', then @protectTerm i@ is the same
-- 'Iteratee', but it will abort if 'stdout' is a terminal device.
protectTerm :: (Nullable s, MonadIO m) => Iteratee s m a -> Iteratee s m a
protectTerm itr = do
    t <- liftIO $ hIsTerminalDevice stdout
    if t then err else itr
  where
    err = error "cowardly refusing to write binary data to terminal"

-- | A simple progress indicator that prints the number of records.
progressNum :: (MonadIO m, Nullable s, NullPoint s, ListLike s a)
            => String -> (String -> IO ()) -> Enumeratee s s m b
progressNum msg put = eneeCheckIfDonePass (icont . go 0)
  where
    go !_ k (EOF   mx) = idone (liftI k) (EOF mx)
    go !n k (Chunk as) = do let !n' = n + LL.length as
                            when (n `div` 65536 /= n' `div` 65536) . liftIO .
                                    put $ "\27[K" ++ msg ++ showNum n ++ "\r"
                            eneeCheckIfDonePass (icont . go n') . k $ Chunk as


-- A very simple queue data type.
-- Invariants: q = QQ l f b --> l == length f + length b
--                          --> l == 0 || not (null f)

data QQ a = QQ !Int [a] [a]

emptyQ :: QQ a
emptyQ = QQ 0 [] []

lengthQ :: QQ a -> Int
lengthQ (QQ l _ _) = l

pushQ :: a -> QQ a -> QQ a
pushQ a (QQ l [] b) = QQ (l+1) (reverse (a:b)) []
pushQ a (QQ l  f b) = QQ (l+1) f (a:b)

popQ :: QQ a -> Maybe (a, QQ a)
popQ (QQ l (a:[]) b) = Just (a, QQ (l-1) (reverse b) [])
popQ (QQ l (a:fs) b) = Just (a, QQ (l-1) fs b)
popQ (QQ _ [    ] _) = Nothing

cancelAll :: MonadIO m => QQ (Async a) -> m ()
cancelAll (QQ _ ff bb) = liftIO $ mapM_ cancel (ff ++ bb)

data ParseError = ParseError {errorContexts :: [String], errorMessage :: String}
    deriving (Show, Typeable)

instance Exception ParseError

-- | A function to convert attoparsec 'Parser's into 'Iteratee's.
parserToIteratee :: (Monad m) => A.Parser a -> Iteratee S.ByteString m a
parserToIteratee p = icont (f (A.parse p)) Nothing
  where
    f k (EOF Nothing) =
        case A.feed (k S.empty) S.empty of
          A.Fail _ err dsc            -> throwErr (toException $ ParseError err dsc)
          A.Partial _                 -> throwErr (toException EofException)
          A.Done rest v | S.null rest -> idone v (EOF Nothing)
                           | otherwise   -> idone v (Chunk rest)
    f _ (EOF (Just e)) = throwErr e
    f k (Chunk s)
        | S.null s = icont (f k) Nothing
        | otherwise =
            case k s of
              A.Fail _ err dsc -> throwErr (toException $ ParseError err dsc)
              A.Partial k'     -> icont (f k') Nothing
              A.Done rest v    -> idone v (Chunk rest)


-- | Equivalent to @joinI $ takeStream n $ stream2vector@, but more
-- efficient.
stream2vectorN :: (MonadIO m, ListLike s a, Nullable s, VG.Vector v a) => Int -> Iteratee s m (v a)
stream2vectorN n = do
    mv <- liftIO $ VM.new n
    l <- go mv 0
    liftIO $ VG.unsafeFreeze $ VM.take l mv
  where
    go mv i
        | i == n    = return n
        | otherwise =
            I.tryHead >>= \x -> case x of
                Nothing -> return i
                Just  a -> liftIO (VM.write mv i a) >> go mv (i+1)

-- | Reads the whole stream into a 'VG.Vector'.
stream2vector :: (MonadIO m, ListLike s a, Nullable s, VG.Vector v a) => Iteratee s m (v a)
stream2vector = liftIO (VM.new 1024) >>= go 0
  where
    go !i !mv = I.tryHead >>= \x -> case x of
                  Nothing -> liftIO $ VG.unsafeFreeze $ VM.take i mv
                  Just  a -> do mv' <- if VM.length mv == i then liftIO (VM.grow mv (VM.length mv)) else return mv
                                when (i `rem` 0x10000 == 0) $ liftIO performGC
                                liftIO $ VM.write mv' i a
                                go (i+1) mv'

withFileFd :: (MonadIO m, MonadMask m) => FilePath -> (Fd -> m a) -> m a
withFileFd filepath iter = bracket
    (liftIO $ openFd filepath ReadOnly Nothing defaultFileFlags)
    (liftIO . closeFd) iter

