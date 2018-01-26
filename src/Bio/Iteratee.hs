module Bio.Iteratee (
    iGetString,
    iterLoop,
    iLookAhead,

    protectTerm,
    parMapChunksIO,
    parRunIO,
    progressGen,
    progressNum,
    progressPos,

    ($==),
    MonadIO, MonadMask,
    lift, liftIO,
    stdin, stdout, stderr,

    enumAuxFile,
    enumInputs,
    enumDefaultInputs,

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

    module Bio.Iteratee.Bytes,
    module Bio.Iteratee.IO,
    module Bio.Iteratee.Iteratee,
    module Bio.Iteratee.List
        ) where

import Bio.Bam.Header
import Bio.Iteratee.Base
import Bio.Iteratee.Bytes
import Bio.Iteratee.IO
import Bio.Iteratee.Iteratee
import Bio.Iteratee.List
import Bio.Prelude
import Bio.Util.Numeric                     ( showNum )
import Control.Concurrent.Async             ( Async, async, wait, cancel )
import Control.Monad.Catch                  ( MonadMask(..) )
import Control.Monad.IO.Class
import Control.Monad.Trans.Class
import System.IO                            ( hIsTerminalDevice )

import qualified Control.Monad.Catch            as CMC
import qualified Data.Attoparsec.ByteString     as A
import qualified Data.ByteString.Char8          as S
import qualified Data.Vector.Generic            as VG
import qualified Data.Vector.Generic.Mutable    as VM

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
iGetString :: Int -> Iteratee S.ByteString m S.ByteString
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
iterLoop it a = do e <- isFinished
                   if e then return a
                        else it a >>= iterLoop it
infixl 1 $==
{-# INLINE ($==) #-}
-- | Compose an 'Enumerator'' with an 'Enumeratee', giving a new
-- 'Enumerator''.
($==) :: Monad m => Enumerator' hdr input m (Iteratee output m result)
                 -> Enumeratee      input             output m result
                 -> Enumerator' hdr                   output m result
($==) enum enee iter = run =<< enum (enee . iter)

-- | Merge two 'Enumerator''s into one.  The header provided by the
-- inner 'Enumerator'' is passed to the output iterator, the header
-- provided by the outer 'Enumerator'' is passed to the merging iteratee

{-# INLINE mergeEnums' #-}
mergeEnums' :: (Nullable s2, Nullable s1, Monad m)
            => Enumerator' hi s1 m a                            -- ^ inner enumerator
            -> Enumerator' ho s2 (Iteratee s1 m) a              -- ^ outer enumerator
            -> (ho -> Enumeratee  s2 s1 (Iteratee s1 m) a)      -- ^ merging enumeratee
            -> Enumerator' hi s1 m a
mergeEnums' e1 e2 etee i = e1 $ \hi -> e2 (\ho -> joinI . etee ho $ ilift lift (i hi)) >>= run

type Enumerator' h eo m b = (h -> Iteratee eo m b) -> m (Iteratee eo m b)
type Enumeratee' h ei eo m b = (h -> Iteratee eo m b) -> Iteratee ei m (Iteratee eo m b)

enumAuxFile :: (MonadIO m, MonadMask m) => FilePath -> Iteratee S.ByteString m a -> m a
enumAuxFile fp it = liftIO (findAuxFile fp) >>= \f -> enumFile defaultBufSize f it >>= run

enumDefaultInputs :: (MonadIO m, MonadMask m) => Enumerator S.ByteString m a
enumDefaultInputs it0 = liftIO getArgs >>= flip enumInputs it0

enumInputs :: (MonadIO m, MonadMask m) => [FilePath] -> Enumerator S.ByteString m a
enumInputs [] = enumFd defaultBufSize stdInput
enumInputs xs = go xs
  where go ("-":fs) = enumFd defaultBufSize stdInput >=> go fs
        go ( f :fs) = enumFile defaultBufSize f >=> go fs
        go [      ] = return

data Ordering' a = Less | Equal a | NotLess

mergeSortStreams :: Monad m => (a -> a -> Ordering' a) -> Enumeratee [a] [a] (Iteratee [a] m) b
mergeSortStreams comp = eneeCheckIfDone step
  where
    step out = peekStream >>= \mx -> lift peekStream >>= \my -> case (mx, my) of
        (Just x, Just y) -> case x `comp` y of
            Less    -> do dropStream 1 ;                       eneeCheckIfDone step . out $ Chunk [x]
            NotLess -> do                lift (dropStream 1) ; eneeCheckIfDone step . out $ Chunk [y]
            Equal z -> do dropStream 1 ; lift (dropStream 1) ; eneeCheckIfDone step . out $ Chunk [z]

        (Just  x, Nothing) -> do       dropStream 1  ; eneeCheckIfDone step . out $ Chunk [x]
        (Nothing, Just  y) -> do lift (dropStream 1) ; eneeCheckIfDone step . out $ Chunk [y]
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
    go' !qq k (EOF  mx) = do a <- liftIO (async (f emptyP))
                             goE mx (pushQ a qq) k Nothing
    go' !qq k (Chunk c) = do a <- liftIO (async (f c))
                             go (pushQ a qq) k Nothing

    -- input ended, empty the queue
    goE  _ !qq k (Just e) = cancelAll qq >> icont (go' emptyQ k) (Just e)
    goE mx !qq k Nothing = case popQ qq of
        Nothing      -> idone (liftI k) (EOF mx)
        Just (a,qq') -> liftIO (wait a) >>= eneeCheckIfDonePass (goE mx qq') . k . Chunk

parRunIO :: MonadIO m => Int -> Enumeratee [IO a] a m b
parRunIO np = eneeCheckIfDonePass (go emptyQ)
  where
    -- check if the queue is full
    go !qq k (Just  e) = cancelAll qq >> icont (go' emptyQ k) (Just e)
    go !qq k  Nothing  = case popQ qq of
        Just (a,qq') | lengthQ qq == np -> liftIO (wait a) >>= eneeCheckIfDonePass (go qq') . k . Chunk
        _                               -> liftI $ go' qq k

    -- we have room for input
    go' !qq k (Chunk (c:cs)) = liftIO (async c) >>= \a -> go' (pushQ a qq) k (Chunk cs)
    go' !qq k (Chunk [    ]) = go qq k Nothing
    go' !qq k (EOF       mx) = goE mx qq k Nothing

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

-- | A general progress indicator that prints some message after a set
-- number of records have passed through.
progressGen :: MonadIO m
            => (Int -> a -> String) -> Int -> (String -> IO ()) -> Enumeratee [a] [a] m b
progressGen msg sz put = eneeCheckIfDonePass (icont . go 0)
  where
    go !_ k (EOF   mx) = idone (liftI k) (EOF mx)
    go !n k (Chunk as)
        | null as    = liftI $ go n k
        | otherwise  = let !n' = n + length as
                       in when (n' `div` sz /= n `div` sz) (liftIO . put $
                                "\27[K" ++ msg n' (head as) ++ "\r")
                          `ioBind_` eneeCheckIfDonePass (icont . go n') (k $ Chunk as)

-- | A simple progress indicator that prints the number of records.
progressNum :: MonadIO m
            => String -> Int -> (String -> IO ()) -> Enumeratee [a] [a] m b
progressNum msg = progressGen (\n _ -> msg ++ " " ++ showNum n)

-- | A simple progress indicator that prints a position every set number
-- of passed records.
progressPos :: MonadIO m
            => (a -> (Refseq, Int)) -> String -> Refs
            -> Int -> (String -> IO ()) -> Enumeratee [a] [a] m b
progressPos f msg refs =
    progressGen $ \_ a -> let (!rs1, !po1) = f a
                              !nm = unpack . sq_name $ getRef refs rs1
                          in msg ++ " " ++ nm ++ ":" ++ showNum po1

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
popQ (QQ _ [    ] _) = Nothing
popQ (QQ l [ a  ] b) = Just (a, QQ (l-1) (reverse b) [])
popQ (QQ l (a:fs) b) = Just (a, QQ (l-1) fs b)

cancelAll :: MonadIO m => QQ (Async a) -> m ()
cancelAll (QQ _ ff bb) = liftIO $ mapM_ cancel (ff ++ bb)

data ParseError = ParseError {errorContexts :: [String], errorMessage :: String}
    deriving (Show, Typeable)

instance Exception ParseError

-- | A function to convert attoparsec 'Parser's into 'Iteratee's.
parserToIteratee :: A.Parser a -> Iteratee S.ByteString m a
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
stream2vectorN :: (MonadIO m, VG.Vector v a) => Int -> Iteratee [a] m (v a)
stream2vectorN n = do
    mv <- liftIO $ VM.new n
    l <- go mv 0
    liftIO $ VG.unsafeFreeze $ VM.take l mv
  where
    go mv i
        | i == n    = return n
        | otherwise =
            tryHead >>= \case
                Nothing -> return i
                Just  a -> liftIO (VM.write mv i a) >> go mv (i+1)

-- | Reads the whole stream into a 'VG.Vector'.
stream2vector :: (MonadIO m, VG.Vector v a) => Iteratee [a] m (v a)
stream2vector = liftIO (VM.new 1024) >>= go 0
  where
    go !i !mv = tryHead >>= \case
                  Nothing -> liftIO $ VG.unsafeFreeze $ VM.take i mv
                  Just  a -> do mv' <- if VM.length mv == i then liftIO (VM.grow mv (VM.length mv)) else return mv
                                when (i `rem` 0x10000 == 0) $ liftIO performGC
                                liftIO $ VM.write mv' i a
                                go (i+1) mv'

withFileFd :: (MonadIO m, MonadMask m) => FilePath -> (Fd -> m a) -> m a
withFileFd filepath = CMC.bracket
    (liftIO $ openFd filepath ReadOnly Nothing defaultFileFlags)
    (liftIO . closeFd)

