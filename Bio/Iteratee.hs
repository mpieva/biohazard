{-# LANGUAGE PatternGuards, BangPatterns #-}

-- | Basically a reexport of @Data.Iteratee@ less the names that clash
-- with @Prelude@ plus a handful of utilities.
module Bio.Iteratee (
    i'groupBy,
    i'groupOn,
    i'getString,
    i'lookAhead,
    i'filterM,
    ($^),
    ListLike,

    module Data.Iteratee.Binary,
    module Data.Iteratee.Char,
    module Data.Iteratee.IO,
    module Data.Iteratee.Iteratee,
    module Data.Iteratee.Parallel
                    ) where

import Control.Monad
import Control.Monad.Trans.Class
import Data.Iteratee.Binary
import Data.Iteratee.Char
import Data.Iteratee.IO
import Data.Iteratee.Iteratee
import Data.Iteratee.Parallel
import Data.Monoid
import Data.ListLike ( ListLike )

import qualified Data.ListLike as LL
import qualified Data.ByteString as S

-- | Grouping on @Iteratee@s.  @i'groupOn proj inner outer@ executes
-- @inner (proj e)@, where @e@ is the first input element, to obtain an
-- @Iteratee i@, then passes elements @e@ to @i@ as long as @proj e@
-- produces the same result.  If @proj e@ changes or the input ends, the
-- pair of @proj e@ and the result of @run i@ is passed to @outer@.  At
-- end of input, the resulting @outer@ is returned.
i'groupOn :: (Monad m, LL.ListLike l e, Eq t1, NullPoint l, Nullable l)
          => (e -> t1)
          -> (t1 -> m (Iteratee l m t2))
          -> Enumeratee l [(t1, t2)] m a
i'groupOn proj inner = eneeCheckIfDone (liftI . step)
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
i'groupBy :: (Monad m, LL.ListLike l t, NullPoint l, Nullable l)
          => (t -> t -> Bool)
          -> m (Iteratee l m t2)
          -> Enumeratee l [t2] m a
i'groupBy cmp inner = eneeCheckIfDone (liftI . step)
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


-- | Run an Iteratee, collect the input.  When it finishes, return the
-- result along with *all* input.  Effectively allows lookahead.  Be
-- careful, this will eat memory if the @Iteratee@ doesn't return
-- speedily.
i'lookAhead :: Monoid s => Iteratee s m a -> Iteratee s m a
i'lookAhead = go mempty
  where 
    go acc it = Iteratee $ \od oc -> runIter it (\x _ -> od x (Chunk acc)) (oc . step acc)
    
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

infixl 2 $^
-- | Compose an @Enumertator@ with an @Enumeratee@, giving a new
-- @Enumerator@.
($^) :: Monad m => Enumerator input m (Iteratee output m result)
                -> Enumeratee input output m result 
                -> Enumerator output m result
($^) enum enee iter = run =<< enum (enee iter)

-- | Apply a monadic filter predicate to an @Iteratee@.
i'filterM :: Monad m => (a -> m Bool) -> Enumeratee [a] [a] m b
i'filterM k = eneeCheckIfDone (liftI . step)
  where
    step it (Chunk xs) = lift (filterM k xs) >>=
                         eneeCheckIfDone (liftI . step) . it . Chunk
    step it stream     = idone (liftI it) stream
 
