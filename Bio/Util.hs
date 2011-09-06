{-# LANGUAGE PatternGuards, BangPatterns #-}
module Bio.Util (
    groupBy,
    groupOn,
    getString,
    lookAheadI,
    R(..)
                ) where

import Control.Monad.Trans.Class
import Data.Iteratee hiding ( groupBy, peek )
import Data.Monoid
import Data.Word ( Word8 )
import Foreign.Ptr ( plusPtr )
import Foreign.Storable ( peek )

import qualified Data.ListLike as LL
import qualified Data.ByteString as S

import Data.ByteString.Unsafe
import Data.ByteString.Internal

-- ^ Stuff that didn't quite fit anywhere else.

-- | Grouping on @Iteratee@s.  @groupOn proj inner outer@ executes
-- @inner (proj e)@, where @e@ is the first input element, to obtain an
-- @Iteratee i@, then passes elements @e@ to @i@ as long as @proj e@
-- produces the same result.  If @proj e@ changes or the input ends, the
-- pair of @proj e@ and the result of @run i@ is passed to @outer@.  At
-- end of input, the resulting @outer@ is returned.

groupOn :: (Monad m, LL.ListLike l e, Eq t1, NullPoint l, Nullable l)
        => (e -> t1)
        -> (t1 -> m (Iteratee l m t2))
        -> Enumeratee l [(t1, t2)] m a
groupOn proj inner = liftI . step
  where
    step outer   (EOF   mx) = idone outer (EOF mx)
    step outer c@(Chunk as)
        | LL.null as = liftI $ step outer
        | otherwise  = lift (inner x) >>= \i -> step' x i outer c
            where !x = proj (LL.head as)

    step' c it outer (Chunk as)
        | LL.null as = liftI $ step' c it outer
        | (l,r) <- LL.span ((==) c . proj) as, not (LL.null l) =
            lift (enumPure1Chunk l it) >>= \it' -> step' c it' outer (Chunk r)

    step' c it outer str = 
        lift (run it >>= \b -> enumPure1Chunk [(c,b)] outer) >>=
        \outer' -> step outer' str


-- | Grouping on @Iteratee@s.  @groupBy cmp inner outer@ executes
-- @inner@ to obtain an @Iteratee i@, then passes elements @e@ to @i@ as
-- long as @cmp e' e@, where @e'@ is some preceeding element, is true.
-- Else, the result of @run i@ is passed to @outer@ and @groupBy@
-- restarts.  At end of input, the resulting @outer@ is returned.

groupBy :: (Monad m, LL.ListLike l t, NullPoint l, Nullable l)
        => (t -> t -> Bool)
        -> m (Iteratee l m t2)
        -> Enumeratee l [t2] m a
groupBy cmp inner = liftI . step
  where
    step outer    (EOF   mx) = idone outer (EOF mx)
    step outer  c@(Chunk as)
        | LL.null as = liftI $ step outer
        | otherwise  = lift inner >>= \i -> step' (LL.head as) i outer c

    step' c it outer (Chunk as)
        | LL.null as = liftI $ step' c it outer
        | (l,r) <- LL.span (cmp c) as, not (LL.null l) =
            lift (enumPure1Chunk l it) >>= \it' -> step' (LL.head l) it' outer (Chunk r)

    step' _ it outer str = 
        lift (run it >>= \b -> enumPure1Chunk [b] outer) >>=
        \outer' -> step outer' str


-- | Run an Iteratee, collect the input.  When it finishes, return the
-- result along with *all* input.  Effectively allows lookahead.  Be
-- careful, this will eat memory if the @Iteratee@ doesn't return
-- speedily.
lookAheadI :: Monoid s => Iteratee s m a -> Iteratee s m a
lookAheadI = go mempty
  where 
    go acc it = Iteratee $ \od oc -> runIter it (\x _ -> od x (Chunk acc)) (oc . step acc)
    
    step acc k c@(Chunk str) = go (acc `mappend` str) (k c)
    step acc k c@(EOF     _) = Iteratee $ \od1 -> runIter (k c) (\x _ -> od1 x (Chunk acc))
                                      

-- | Collects a string of a given length.  Don't use this for long
-- strings, use @Data.Iteratee.ListLike.take@ instead.
getString :: Monad m => Int -> Iteratee S.ByteString m S.ByteString
getString 0 = idone S.empty (Chunk S.empty)
getString n = liftI $ step [] 0
  where
    step acc l c@(EOF _) = icont (step acc l) (Just $ setEOF c)
    step acc l (Chunk c) | l + S.length c >= n = let r = S.concat . reverse $ S.take (n-l) c : acc
                                                 in idone r (Chunk $ S.drop (n-l) c)
                         | otherwise           = liftI $ step (c:acc) (l + S.length c)


-- Occasionally faster method to compare @Seqid@s, by starting at the
-- end.  This makes sense because people tend to name their reference
-- sequences like "contig_xxx", so comparing the beginning isn't really
-- helpful.  Same goes for query sequences, which tend to start will
-- longish runnames.
newtype R = R S.ByteString

instance Ord R where compare = compare_R
instance Eq  R where a == b = case compare_R a b of EQ -> True ; _ -> False

compare_R :: R -> R -> Ordering
compare_R (R a) (R b) = inlinePerformIO $
                        unsafeUseAsCStringLen a $ \(pa,la) -> 
                        unsafeUseAsCStringLen b $ \(pb,lb) ->
                        case compare la lb of LT -> return LT
                                              GT -> return GT
                                              EQ -> go (pa `plusPtr` (la-1)) (pb `plusPtr` (lb-1)) la
    where
        go !_ !_ 0 = return EQ
        go  p  q n = do x <- peek p :: IO Word8
                        y <- peek q :: IO Word8
                        case compare x y of
                              LT -> return LT
                              GT -> return GT
                              EQ -> go (p `plusPtr` 1) (q `plusPtr` 1) (n-1)

