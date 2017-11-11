{-# LANGUAGE TupleSections, ScopedTypeVariables #-}

-- |Monadic Iteratees:
-- incremental input parsers, processors and transformers
--
-- This module provides many basic iteratees from which more complicated
-- iteratees can be built.  In general these iteratees parallel those in
-- @Data.List@, with some additions.

module Bio.Iteratee.List (
  -- * Iteratees
  -- ** Iteratee Utilities
  isFinished
  ,stream2list
  ,stream2stream
  -- ** Basic Iteratees
  ,dropWhileStream
  ,dropStream
  ,headStream
  ,tryHead
  ,lastStream
  ,heads
  ,peekStream
  ,roll
  ,lengthStream
  ,chunkLength
  ,takeFromChunk
  -- ** Nested iteratee combinators
  ,breakStream
  ,breakE
  ,takeStream
  ,takeUpTo
  ,takeWhileE
  ,mapStream
  ,concatMapStream
  ,concatMapStreamM
  ,mapMaybeStream
  ,filterStream
  ,filterStreamM
  ,groupStreamBy
  ,groupStreamOn
  ,mergeStreams
  ,mergeByChunks
  -- ** Folds
  ,foldStream
  -- * Enumerators
  -- ** Basic enumerators
  ,enumPureNChunk
  -- ** Enumerator Combinators
  ,enumWith
  ,zipStreams
  ,zipStreams3
  ,zipStreams4
  ,zipStreams5
  ,sequenceStreams_
  ,countConsumed
  -- ** Monadic functions
  ,mapStreamM
  ,mapStreamM_
  ,foldStreamM
  -- * Re-exported modules
  ,module Bio.Iteratee.Iteratee
)
where

import Bio.Iteratee.Iteratee
import Bio.Prelude
import Control.Monad.Trans.Class

-- import qualified Data.ByteString          as B


-- Useful combinators for implementing iteratees and enumerators

-- | Check if a stream has received 'EOF'.
isFinished :: Nullable s => Iteratee s m Bool
isFinished = liftI check
  where
  check c@(Chunk xs)
    | nullC xs    = liftI check
    | otherwise   = idone False c
  check s@(EOF _) = idone True s
{-# INLINE isFinished #-}

-- ------------------------------------------------------------------------
-- Primitive iteratees

-- |Read a stream to the end and return all of its elements as a list.
-- This iteratee returns all data from the stream *strictly*.
stream2list :: Monad m => Iteratee [el] m [el]
stream2list = liftM concat getChunks
{-# INLINE stream2list #-}

-- |Read a stream to the end and return all of its elements as a stream.
-- This iteratee returns all data from the stream *strictly*.
stream2stream :: (Monad m, Nullable s, Monoid s) => Iteratee s m s
stream2stream = liftM mconcat getChunks
{-# INLINE stream2stream #-}


-- ------------------------------------------------------------------------
-- Parser combinators

-- |Attempt to read the next element of the stream and return it
-- Raise a (recoverable) error if the stream is terminated.
--
-- The analogue of @List.head@
--
-- Because @head@ can raise an error, it shouldn't be used when constructing
-- iteratees for @convStream@.  Use @tryHead@ instead.
headStream :: Iteratee [el] m el
headStream = liftI step
  where
  step (Chunk [     ]) = icont step Nothing
  step (Chunk (hd:tl)) = idone hd (Chunk tl)
  step stream          = icont step (Just (setEOF stream))
{-# INLINE headStream #-}

-- | Similar to @headStream@, except it returns @Nothing@ if the stream
-- is terminated.
tryHead :: Iteratee [el] m (Maybe el)
tryHead = liftI step
  where
  step (Chunk [     ]) = liftI step
  step (Chunk (hd:tl)) = idone (Just hd) (Chunk tl)
  step stream          = idone Nothing stream
{-# INLINE tryHead #-}

-- |Attempt to read the last element of the stream and return it
-- Raise a (recoverable) error if the stream is terminated
--
-- The analogue of @List.last@
lastStream :: Iteratee [el] m el
lastStream = liftI (step Nothing)
  where
  step l (Chunk xs)
    | nullC xs     = liftI (step l)
    | otherwise    = liftI $ step (Just $ last xs)
  step l s@(EOF _) = case l of
    Nothing -> icont (step l) . Just . setEOF $ s
    Just x  -> idone x s
{-# INLINE lastStream #-}


-- |Given a sequence of characters, attempt to match them against
-- the characters on the stream.  Return the count of how many
-- characters matched.  The matched characters are removed from the
-- stream.
-- For example, if the stream contains 'abd', then (heads 'abc')
-- will remove the characters 'ab' and return 2.
heads :: (Monad m, Eq el) => [el] -> Iteratee [el] m Int
heads st | nullC st = return 0
heads st = loopE 0 st
  where
  loopE cnt xs
    | nullC xs  = return cnt
    | otherwise = liftI (step cnt xs)
  step cnt str (Chunk [])          = liftI (step cnt str)
  step cnt [ ] stream              = idone cnt stream
  step cnt (y:ys) s@(Chunk (x:xs))
    | y == x    = step (succ cnt) ys (Chunk xs)
    | otherwise = idone cnt s
  step cnt _ stream         = idone cnt stream
{-# INLINE heads #-}


-- |Look ahead at the next element of the stream, without removing
-- it from the stream.
-- Return @Just c@ if successful, return @Nothing@ if the stream is
-- terminated by 'EOF'.
peekStream :: Iteratee [el] m (Maybe el)
peekStream = liftI step
  where
    step   (Chunk [   ]) = liftI step
    step s@(Chunk (x:_)) = idone (Just x) s
    step stream          = idone Nothing stream
{-# INLINE peekStream #-}

-- | Return a chunk of @t@ elements length while consuming @d@ elements
--   from the stream.  Useful for creating a 'rolling average' with
--  'convStream'.
roll
  :: Monad m
  => Int  -- ^ length of chunk (t)
  -> Int  -- ^ amount to consume (d)
  -> Iteratee [el] m [[el]]
roll t d | t > d  = liftI step
  where
    step (Chunk vec)
      | length vec >= t =
          idone [take t vec] (Chunk $ drop d vec)
      | null vec        = liftI step
      | otherwise          = liftI (step' vec)
    step stream            = idone empty stream
    step' v1 (Chunk vec)   = step . Chunk $ v1 `mappend` vec
    step' v1 stream        = idone [v1] stream
roll t d = do r <- joinI (takeStream t stream2stream)
              dropStream (d-t)
              return [r]
  -- d is >= t, so this version works
{-# INLINE roll #-}


-- |Drop n elements of the stream, if there are that many.
--
-- The analogue of @List.drop@
dropStream :: Int -> Iteratee [el] m ()
dropStream 0  = idone () (Chunk emptyP)
dropStream n' = liftI (step n')
  where
    step n (Chunk str)
      | length str < n = liftI (step (n - length str))
      | otherwise         = idone () (Chunk (drop n str))
    step _ stream         = idone () stream
{-# INLINE dropStream #-}

-- |Skip all elements while the predicate is true.
--
-- The analogue of @List.dropWhile@
dropWhileStream :: (el -> Bool) -> Iteratee [el] m ()
dropWhileStream p = liftI step
  where
    step (Chunk str)
      | null rest = liftI step
      | otherwise    = idone () (Chunk rest)
      where
        rest = dropWhile p str
    step stream      = idone () stream
{-# INLINE dropWhileStream #-}


-- | Return the total length of the remaining part of the stream.
--
-- This forces evaluation of the entire stream.
--
-- The analogue of @List.length@
lengthStream :: Num a => Iteratee [el] m a
lengthStream = liftI (step 0)
  where
    step !i (Chunk xs) = liftI (step $ i + fromIntegral (length xs))
    step !i stream     = idone i stream
{-# INLINE lengthStream #-}

-- | Get the length of the current chunk, or @Nothing@ if 'EOF'.
--
-- This function consumes no input.
chunkLength :: Iteratee [el] m (Maybe Int)
chunkLength = liftI step
 where
  step s@(Chunk xs) = idone (Just $ length xs) s
  step stream       = idone Nothing stream
{-# INLINE chunkLength #-}

-- | Take @n@ elements from the current chunk, or the whole chunk if
-- @n@ is greater.
takeFromChunk :: Int -> Iteratee [el] m [el]
takeFromChunk n | n <= 0 = idone emptyP (Chunk emptyP)
takeFromChunk n = liftI step
 where
  step (Chunk xs) = let (h,t) = splitAt n xs in idone h $ Chunk t
  step stream     = idone emptyP stream
{-# INLINE takeFromChunk #-}

-- |Takes an element predicate and returns the (possibly empty) prefix of
-- the stream.  None of the characters in the string satisfy the character
-- predicate.
-- If the stream is not terminated, the first character of the remaining stream
-- satisfies the predicate.
--
-- N.B. 'breakE' should be used in preference to @breakStream@.
-- @breakStream@ will retain all data until the predicate is met, which may
-- result in a space leak.
--
-- The analogue of @List.break@

breakStream :: (el -> Bool) -> Iteratee [el] m [el]
breakStream cpred = icont (step mempty) Nothing
  where
    step bfr (Chunk str)
      | null str          =  icont (step bfr) Nothing
      | otherwise         =  case break cpred str of
        (str', tail')
          | null tail'    -> icont (step (bfr `mappend` str)) Nothing
          | otherwise     -> idone (bfr `mappend` str') (Chunk tail')
    step bfr stream       =  idone bfr stream
{-# INLINE breakStream #-}

-- ---------------------------------------------------
-- The converters show a different way of composing two iteratees:
-- `vertical' rather than `horizontal'

-- |Takes an element predicate and an iteratee, running the iteratee
-- on all elements of the stream until the predicate is met.
--
-- the following rule relates @break@ to @breakE@
-- @break@ pred === @joinI@ (@breakE@ pred stream2stream)
--
-- @breakE@ should be used in preference to @break@ whenever possible.
breakE :: (el -> Bool) -> Enumeratee [el] [el] m a
breakE cpred = eneeCheckIfDonePass (icont . step)
 where
  step k (Chunk s)
      | null s  = liftI (step k)
      | otherwise  = case break cpred s of
        (str', tail')
          | null tail'    -> eneeCheckIfDonePass (icont . step) . k $ Chunk str'
          | otherwise     -> idone (k $ Chunk str') (Chunk tail')
  step k stream           =  idone (liftI k) stream
{-# INLINE breakE #-}

-- |Read n elements from a stream and apply the given iteratee to the
-- stream of the read elements. Unless the stream is terminated early, we
-- read exactly n elements, even if the iteratee has accepted fewer.
--
-- The analogue of @List.take@
takeStream ::
  Monad m
  => Int   -- ^ number of elements to consume
  -> Enumeratee [el] [el] m a
takeStream n' iter
 | n' <= 0   = return iter
 | otherwise = Iteratee $ \od oc -> runIter iter (on_done od oc) (on_cont od oc)
  where
    on_done od oc x _ = runIter (dropStream n' >> return (return x)) od oc
    on_cont od oc k Nothing = if n' == 0 then od (liftI k) (Chunk mempty)
                                 else runIter (liftI (step n' k)) od oc
    on_cont od oc _ (Just e) = runIter (dropStream n' >> throwErr e) od oc
    step n k (Chunk str)
      | null str           = liftI (step n k)
      | length str <= n    = takeStream (n - length str) $ k (Chunk str)
      | otherwise          = idone (k (Chunk s1)) (Chunk s2)
      where (s1, s2) = splitAt n str
    step _n k stream       = idone (liftI k) stream
{-# INLINE takeStream #-}

-- |Read n elements from a stream and apply the given iteratee to the
-- stream of the read elements. If the given iteratee accepted fewer
-- elements, we stop.
-- This is the variation of 'takeStream' with the early termination
-- of processing of the outer stream once the processing of the inner stream
-- finished early.
--
-- Iteratees composed with 'takeUpTo' will consume only enough elements to
-- reach a done state.  Any remaining data will be available in the outer
-- stream.
--
-- > > let iter = do
-- > h <- joinI $ takeUpTo 5 I.head
-- > t <- stream2list
-- > return (h,t)
-- >
-- > > enumPureNChunk [1..10::Int] 3 iter >>= run >>= print
-- > (1,[2,3,4,5,6,7,8,9,10])
-- >
-- > > enumPureNChunk [1..10::Int] 7 iter >>= run >>= print
-- > (1,[2,3,4,5,6,7,8,9,10])
--
-- in each case, @I.head@ consumes only one element, returning the remaining
-- 4 elements to the outer stream
takeUpTo :: Monad m => Int -> Enumeratee [el] [el] m a
takeUpTo i iter
 | i <= 0    = idone iter (Chunk emptyP)
 | otherwise = Iteratee $ \od oc ->
    runIter iter (onDone od oc) (onCont od oc)
  where
    onDone od oc x str      = runIter (idone (return x) str) od oc
    onCont od oc k Nothing  = if i == 0 then od (liftI k) (Chunk mempty)
                                 else runIter (liftI (step i k)) od oc
    onCont od oc _ (Just e) = runIter (throwErr e) od oc
    step n k (Chunk str)
      | null str       = liftI (step n k)
      | length str < n = takeUpTo (n - length str) $ k (Chunk str)
      | otherwise      =
         -- check to see if the inner iteratee has completed, and if so,
         -- grab any remaining stream to put it in the outer iteratee.
         -- the outer iteratee is always complete at this stage, although
         -- the inner may not be.
         let (s1, s2) = splitAt n str
         in Iteratee $ \od' _ -> do
              res <- runIter (k (Chunk s1)) (\a s  -> return $ Left  (a, s))
                                            (\k' e -> return $ Right (k',e))
              case res of
                Left (a,Chunk s1') -> od' (return a)
                                          (Chunk $ s1' ++ s2)
                Left  (a,s')       -> od' (idone a s') (Chunk s2)
                Right (k',e)       -> od' (icont k' e) (Chunk s2)
    step _ k stream       = idone (liftI k) stream
{-# INLINE takeUpTo #-}


-- |Takes an element predicate and an iteratee, running the iteratee
-- on all elements of the stream while the predicate is met.
--
-- This is preferred to @takeWhile@.
takeWhileE :: (el -> Bool) -> Enumeratee [el] [el] m a
takeWhileE = breakE . (not .)
{-# INLINEABLE takeWhileE #-}

-- | Map a function over an 'Iteratee'.
-- This one is reimplemented and differs from the the one in
-- "Data.Iteratee.ListLike" in so far that it doesn't pass on an 'EOF'
-- received in the input, which is the expected behavior.
mapStream :: (el -> el') -> Enumeratee [el] [el'] m a
mapStream = mapChunks . map
{-# INLINE mapStream #-}

-- | Apply a function to the elements of a stream, concatenate the
-- results into a stream.  No giant intermediate list is produced.
concatMapStream :: Monoid t => (a -> t) -> Enumeratee [a] t m r
concatMapStream = mapChunks . foldMap
{-# INLINE concatMapStream #-}

-- | Apply a monadic function to the elements of a stream, concatenate
-- the results into a stream.  No giant intermediate list is produced.
concatMapStreamM :: Monad m => (a -> m t) -> Enumeratee [a] t m r
concatMapStreamM f = eneeCheckIfDone (liftI . go)
  where
    go k (EOF   mx)              = idone (liftI k) (EOF mx)
    go k (Chunk xs) | null xs    = liftI (go k)
                    | otherwise  = f (head xs) `mBind`
                                   eneeCheckIfDone (flip go (Chunk (tail xs))) . k . Chunk
{-# INLINE concatMapStreamM #-}

mapMaybeStream :: (a -> Maybe b) -> Enumeratee [a] [b] m r
mapMaybeStream = mapChunks . mapMaybe
{-# INLINE mapMaybeStream #-}

-- |Creates an 'enumeratee' with only elements from the stream that
-- satisfy the predicate function.  The outer stream is completely consumed.
--
-- The analogue of @List.filter@
filterStream :: (el -> Bool) -> Enumeratee [el] [el] m a
filterStream p = mapChunks (filter p)
{-# INLINE filterStream #-}

-- | Apply a monadic filter predicate to an 'Iteratee'.
filterStreamM :: Monad m => (a -> m Bool) -> Enumeratee [a] [a] m r
filterStreamM k = mapChunksM (go id)
  where
    go acc [   ] = return $! acc empty
    go acc (h:t) = do p <- k h
                      let acc' = if p then (:) h . acc else acc
                      go acc' t
{-# INLINE filterStreamM #-}

-- | Grouping on 'Iteratee's.  @groupStreamOn proj inner outer@ executes
-- @inner (proj e)@, where @e@ is the first input element, to obtain an
-- 'Iteratee' @i@, then passes elements @e@ to @i@ as long as @proj e@
-- produces the same result.  If @proj e@ changes or the input ends, the
-- pair of @proj e@ and the result of @run i@ is passed to @outer@.  At
-- end of input, the resulting @outer@ is returned.
groupStreamOn :: (Monad m, Eq t1)
              => (e -> t1)
              -> (t1 -> m (Iteratee [e] m t2))
              -> Enumeratee [e] [(t1, t2)] m a
groupStreamOn proj inner = eneeCheckIfDonePass (icont . step)
  where
    step outer   (EOF      mx) = idone (liftI outer) $ EOF mx
    step outer   (Chunk [   ]) = liftI $ step outer
    step outer c@(Chunk (h:_)) = let x = proj h
                                 in lift (inner x) >>= \i -> step' x i outer c

    -- We want to feed a 'Chunk' to the inner 'Iteratee', which might be
    -- finished.  In that case, we would want to abort, but we cannot,
    -- since the outer iteration is still going on.  So instead we
    -- discard data we would have fed to the inner 'Iteratee'.  (Use of
    -- 'enumPure1Chunk' is not appropriate, it would accumulate the
    -- data, just to have it discarded by the 'run' that eventually
    -- happens.

    step' c it outer (Chunk as)
        | null as = liftI $ step' c it outer
        | (l,r) <- span ((==) c . proj) as, not (null l) =
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
groupStreamBy :: Monad m
              => (t -> t -> Bool)
              -> m (Iteratee [t] m t2)
              -> Enumeratee [t] [t2] m a
groupStreamBy cmp inner = eneeCheckIfDonePass (icont . step)
  where
    step outer   (EOF      mx) = idone (liftI outer) $ EOF mx
    step outer   (Chunk [   ]) = liftI $ step outer
    step outer c@(Chunk (h:_)) = lift inner >>= \i -> step' h i outer c

    step' c it outer (Chunk as)
        | null as = liftI $ step' c it outer
        | (l,r) <- span (cmp c) as, not (null l) =
            let od a    _str = idoneM a $ EOF Nothing
                oc k Nothing = return $ k (Chunk l)
                oc k       m = icontM k m
            in lift (runIter it od oc) >>= \it' -> step' (head l) it' outer (Chunk r)

    step' _ it outer str =
        lift (run it) >>= \b -> eneeCheckIfDone (`step` str) . outer $ Chunk [b]


-- | @mergeStreams@ offers another way to nest iteratees: as a monad stack.
-- This allows for the possibility of interleaving data from multiple
-- streams.
--
-- > -- print each element from a stream of lines.
-- > logger :: (MonadIO m) => Iteratee [ByteString] m ()
-- > logger = mapStreamM_ (liftIO . putStrLn . B.unpack)
-- >
-- > -- combine alternating lines from two sources
-- > -- To see how this was derived, follow the types from
-- > -- 'ileaveLines logger' and work outwards.
-- > run =<< enumFile 10 "file1" (joinI $ enumLinesBS $
-- >           ( enumFile 10 "file2" . joinI . enumLinesBS $ joinI
-- >                 (ileaveLines logger)) >>= run)
-- >
-- > ileaveLines :: (Functor m, Monad m)
-- >   => Enumeratee [ByteString] [ByteString] (Iteratee [ByteString] m)
-- >        [ByteString]
-- > ileaveLines = mergeStreams (\l1 l2 ->
-- >    [B.pack "f1:\n\t" ,l1 ,B.pack "f2:\n\t" ,l2 ]
-- >
-- >
--
mergeStreams :: Monad m => (el1 -> el2 -> b) -> Enumeratee [el2] b (Iteratee [el1] m) a
mergeStreams f = convStream $ liftM2 f (lift headStream) headStream
{-# INLINE mergeStreams #-}

-- | A version of mergeStreams which operates on chunks instead of
-- elements.
--
-- mergeByChunks offers more control than 'mergeStreams'.
-- 'mergeStreams' terminates when the first stream terminates, however
-- mergeByChunks will continue until both streams are exhausted.
--
-- 'mergeByChunks' guarantees that both chunks passed to the merge
-- function will have the same number of elements, although that number
-- may vary between calls.
mergeByChunks ::
  Monad m
  => ([el1] -> [el2] -> c3)  -- ^ merge function
  -> ([el1] -> c3)
  -> ([el2] -> c3)
  -> Enumeratee [el2] c3 (Iteratee [el1] m) a
mergeByChunks f f1 f2 = unfoldConvStream iter (0 :: Int)
 where
  iter 1 = (\x -> (1,f1 x)) `liftM` lift getChunk
  iter 2 = (\x -> (2,f2 x)) `liftM` getChunk
  iter _ = do
    ml1 <- lift chunkLength
    ml2 <- chunkLength
    case (ml1, ml2) of
      (Just l1, Just l2) -> do
        let tval = min l1 l2
        c1 <- lift $ takeFromChunk tval
        c2 <- takeFromChunk tval
        return (0, f c1 c2)
      (Just _, Nothing) -> iter 1
      (Nothing, _)      -> iter 2
{-# INLINE mergeByChunks #-}

-- ------------------------------------------------------------------------
-- Folds

-- | Left-associative fold that is strict in the accumulator.
-- This function should be used in preference to 'foldl' whenever possible.
--
-- The analogue of @List.foldl'@.
foldStream :: (a -> el -> a) -> a -> Iteratee [el] m a
foldStream f i = liftI (step i)
  where
    step acc (Chunk xs)
      | null xs = liftI (step acc)
      | otherwise  = liftI (step $! foldl' f acc xs)
    step acc stream = idone acc stream
{-# INLINE foldStream #-}

-- ------------------------------------------------------------------------
-- Zips

-- |Enumerate two iteratees over a single stream simultaneously.
--
-- Compare to @List.zip@.
zipStreams
  :: Monad m
  => Iteratee [el] m a
  -> Iteratee [el] m b
  -> Iteratee [el] m (a, b)
zipStreams x0 y0 = do
    -- need to check if both iteratees are initially finished.  If so,
    -- we don't want to push a chunk which will be dropped
    (a', x') <- lift $ runIter x0 od oc
    (b', y') <- lift $ runIter y0 od oc
    case checkDone a' b' of
      Just (Right (a,b,s))  -> idone (a,b) s  -- 's' may be EOF, needs to stay
      Just (Left (Left a))  -> liftM (a,) y'
      Just (Left (Right b)) -> liftM (,b) x'
      Nothing               -> liftI (step x' y')
  where
    step x y (Chunk xs) | nullC xs = liftI (step x y)
    step x y (Chunk xs) = do
      (a', x') <- lift $ (\i -> runIter i od oc) =<< enumPure1Chunk xs x
      (b', y') <- lift $ (\i -> runIter i od oc) =<< enumPure1Chunk xs y
      case checkDone a' b' of
        Just (Right (a,b,s))  -> idone (a,b) s
        Just (Left (Left a))  -> liftM (a,) y'
        Just (Left (Right b)) -> liftM (,b) x'
        Nothing               -> liftI (step x' y')
    step x y (EOF err) = joinIM $ case err of
      Nothing -> (liftM2.liftM2) (,) (enumEof   x) (enumEof   y)
      Just e  -> (liftM2.liftM2) (,) (enumErr e x) (enumErr e y)

    od a s = return (Just (a, s), idone a s)
    oc k e = return (Nothing    , icont k e)

    checkDone r1 r2 = case (r1, r2) of
      (Just (a, s1), Just (b,s2)) -> Just $ Right (a, b, shorter s1 s2)
      (Just (a, _), Nothing)      -> Just . Left $ Left a
      (Nothing, Just (b, _))      -> Just . Left $ Right b
      (Nothing, Nothing)          -> Nothing

    shorter c1@(Chunk xs) c2@(Chunk ys)
      | length xs < length ys = c1
      | otherwise                   = c2
    shorter e@(EOF _)  _         = e
    shorter _          e@(EOF _) = e
{-# INLINE zipStreams #-}

zipStreams3
  :: Monad m
  => Iteratee [el] m a -> Iteratee [el] m b
  -> Iteratee [el] m c -> Iteratee [el] m (a, b, c)
zipStreams3 a b c = zipStreams a (zipStreams b c) >>=
  \(r1, (r2, r3)) -> return (r1, r2, r3)
{-# INLINE zipStreams3 #-}

zipStreams4
  :: Monad m
  => Iteratee [el] m a -> Iteratee [el] m b
  -> Iteratee [el] m c -> Iteratee [el] m d
  -> Iteratee [el] m (a, b, c, d)
zipStreams4 a b c d = zipStreams a (zipStreams3 b c d) >>=
  \(r1, (r2, r3, r4)) -> return (r1, r2, r3, r4)
{-# INLINE zipStreams4 #-}

zipStreams5
  :: Monad m
  => Iteratee [el] m a -> Iteratee [el] m b
  -> Iteratee [el] m c -> Iteratee [el] m d
  -> Iteratee [el] m e -> Iteratee [el] m (a, b, c, d, e)
zipStreams5 a b c d e = zipStreams a (zipStreams4 b c d e) >>=
  \(r1, (r2, r3, r4, r5)) -> return (r1, r2, r3, r4, r5)
{-# INLINE zipStreams5 #-}

-- | Enumerate over two iteratees in parallel as long as the first iteratee
-- is still consuming input.  The second iteratee will be terminated with EOF
-- when the first iteratee has completed.  An example use is to determine
-- how many elements an iteratee has consumed:
--
-- > snd <$> enumWith (dropWhile (<5)) length
--
-- Compare to @zipStreams@
enumWith
  :: Monad m
  => Iteratee [el] m a
  -> Iteratee [el] m b
  -> Iteratee [el] m (a, b)
enumWith i1 i2 = do
    -- as with zipStreams, first check to see if the initial iteratee is complete,
    -- otherwise data would be dropped.
    -- running the second iteratee as well to prevent a monadic effect mismatch
    -- although I think that would be highly unlikely to happen in common
    -- code
    (a', x') <- lift $ runIter i1 od oc
    (_,  y') <- lift $ runIter i2 od oc
    case a' of
      Just (a, s) -> flip idone s =<< lift (liftM (a,) $ run i2)
      Nothing     -> go x' y'
  where
    od a s = return (Just (a, s), idone a s)
    oc k e = return (Nothing    , icont k e)

    getUsed xs (Chunk ys) = take (length xs - length ys) xs
    getUsed xs (EOF _)    = xs

    go x y = liftI step
      where
        step (Chunk xs) | nullC xs = liftI step
        step (Chunk xs) = do
          (a', x') <- lift $ (\i -> runIter i od oc) =<< enumPure1Chunk xs x
          case a' of
            Just (a, s) -> do
              b <- lift $ run =<< enumPure1Chunk (getUsed xs s) y
              idone (a, b) s
            Nothing        -> lift (enumPure1Chunk xs y) >>= go x'
        step (EOF err) = joinIM $ case err of
          Nothing -> (liftM2.liftM2) (,) (enumEof   x) (enumEof   y)
          Just e  -> (liftM2.liftM2) (,) (enumErr e x) (enumErr e y)
{-# INLINE enumWith #-}

-- |Enumerate a list of iteratees over a single stream simultaneously
-- and discard the results. This is a different behavior than Prelude's
-- sequence_ which runs iteratees in the list one after the other.
--
-- Compare to @Prelude.sequence_@.
sequenceStreams_
  :: Monad m
  => [Iteratee [el] m a]
  -> Iteratee [el] m ()
sequenceStreams_ = self
  where
    self is = liftI step
      where
        step (Chunk xs) | null xs = liftI step
        step s@(Chunk _) = do
          -- give a chunk to each iteratee
          is'  <- lift $ mapM (enumChunk s) is
          -- filter done iteratees
          (done, notDone) <- lift $ partition fst `liftM` mapM enumCheckIfDone is'
          if null notDone
            then idone () <=< remainingStream $ map snd done
            else self $ map snd notDone
        step s@(EOF _) = do
          s' <- remainingStream <=< lift $ mapM (enumChunk s) is
          case s' of
            EOF (Just e) -> throwErr e
            _            -> idone () s'

    -- returns the unconsumed part of the stream; "sequenceStreams_ is" consumes as
    -- much of the stream as the iteratee in is that consumes the most; e.g.
    -- sequenceStreams_ [I.head, I.last] consumes whole stream
    remainingStream :: Monad m => [Iteratee [el] m a] -> Iteratee [el] m (Stream [el])
    remainingStream is = lift $
      return . foldl1 shorter <=< mapM (\i -> runIter i od oc) $ is
      where
        od _ s = return s
        oc _ e = return $ case e of
          Nothing -> mempty
          _       -> EOF e

    -- return the shorter one of two streams; errors are propagated with the
    -- priority given to the "left"
    shorter c1@(Chunk xs) c2@(Chunk ys)
      | length xs < length ys = c1
      | otherwise                   = c2
    shorter (EOF e1 ) (EOF e2 ) = EOF (e1 `mplus` e2)
    shorter e@(EOF _) _         = e
    shorter _         e@(EOF _) = e

-- |Transform an iteratee into one that keeps track of how much data it
-- consumes.
countConsumed :: (Monad m, Integral n) => Iteratee [el] m a -> Iteratee [el] m (a, n)
countConsumed i = go 0 (const i) (Chunk emptyP)
  where
    go !n f str@(EOF _) = (, n) `liftM` f str
    go !n f str@(Chunk c) = Iteratee rI
      where
        newLen = n + fromIntegral (length c)
        rI od oc = runIter (f str) onDone onCont
          where
            onDone a str'@(Chunk c') =
                od (a, newLen - fromIntegral (length c')) str'
            onDone a str'@(EOF _) = od (a, n) str'
            onCont f' = oc (go newLen f')
{-# INLINE countConsumed #-}

-- ------------------------------------------------------------------------
-- Enumerators

-- |The pure n-chunk enumerator
-- It passes a given stream of elements to the iteratee in @n@-sized chunks.
enumPureNChunk :: Monad m => [el] -> Int -> Enumerator [el] m a
enumPureNChunk str n iter
  | null str = return iter
  | n > 0       = enum' str iter
  | otherwise   = error $ "enumPureNChunk called with n==" ++ show n
  where
    enum' str' iter'
      | null str' = return iter'
      | otherwise    = let (s1, s2) = splitAt n str'
                           on_cont k Nothing = enum' s2 . k $ Chunk s1
                           on_cont k e = return $ icont k e
                       in runIter iter' idoneM on_cont
{-# INLINE enumPureNChunk #-}

-- ------------------------------------------------------------------------
-- Monadic functions

-- | Maps a monadic function over the elements of the stream and ignores
-- the result.
mapStreamM_ :: Monad m => (el -> m b) -> Iteratee [el] m ()
mapStreamM_ = mapChunksM_ . mapM_
{-# INLINE mapStreamM_ #-}

-- | Maps a monadic function over an 'Iteratee'.
mapStreamM :: Monad m => (el -> m el') -> Enumeratee [el] [el'] m a
mapStreamM = mapChunksM . mapM
{-# INLINE mapStreamM #-}

-- | Folds a monadic function over an 'Iteratee'.
foldStreamM :: Monad m => (b -> a -> m b) -> b -> Iteratee [a] m b
foldStreamM = foldChunksM . foldM
{-# INLINE foldStreamM #-}
