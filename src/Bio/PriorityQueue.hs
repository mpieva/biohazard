{-# LANGUAGE BangPatterns #-}
module Bio.PriorityQueue (
        Sizeable(..),
        PQ_Conf(..),

        PQ,
        withPQ,
        makePQ,
        deletePQ,
        enqueuePQ,
        dequeuePQ,
        getMinPQ,
        peekMinPQ,
        sizePQ
) where

import Data.Binary
import Data.IORef
import qualified Control.Exception as CE
import Prelude

-- | A Priority Queue that can fall back to external storage.
--
-- Note that such a Priority Queue automatically gives rise to an
-- external sorting algorithm:  enqueue everything, dequeue until empty.
--
-- Whatever is to be stored in this queue needs to be in Binary, because
-- it may need to be moved to external storage on demand.  We also need
-- a way to estimate the memory consumption of an enqueued object.  When
-- constructing the queue, the maximum amount of RAM to consume is set.
-- Note that open input streams use memory for buffering, too.
--
-- Enqueued objects are kept in an in memory heap until the memory
-- consumption becomes too high.  At that point, the whole heap is
-- sorted and dumped to external storage.  If necessary, the file to do
-- so is created and kept open.  The newly created stream is added to a
-- heap so that dequeueing objects amounts to performing a merge sort on
-- multiple external streams.  To conserve on file descriptors, we
-- concatenate multiple streams into a single file, then use pread(2) on
-- that as appropriate.  If too many streams are open (how do we set
-- that limit?), we do exactly that:  merge-sort all streams and the
-- in-memory heap into a single new stream.  One file is created for
-- each generation of streams, so that mergind handles streams of
-- roughly equal length.
--
-- XXX  Truth be told, this queue isn't backed externally, and ignores
--      all limits.  It *is* a Priority Queue, though!
--
-- XXX  May want to add callbacks for significant events (new file,
--      massive merge, deletion of file?)
--
-- XXX  Need to track memory consumption of input buffers.
--
-- XXX  Need a way to decide when too many streams are open.  That point
--      is reached when seeking takes about as much time as reading
--      (which depends on buffer size and system characteristics), so
--      that an additional merge pass becomes economical.
--
-- XXX  These will be useful:
--          unix-bytestring:System.Posix.IO.ByteString.fdPread
--          temporary:System.IO.Temp.openBinaryTempFile
--          lz4:Codec.Compression.LZ4

data PQ_Conf = PQ_Conf {
        max_mb :: Int,          -- ^ memory limit
        temp_path :: FilePath   -- ^ path to temporary files (a directory will be created)
        -- functions to report progress go here
    }

newtype PQ a = PQ (IORef (SkewHeap a, Int))

class Sizeable a where usedBytes :: a -> Int

-- | Creates a priority queue.  Note that the priority queue creates
-- files, which will only be cleaned up if deletePQ is called.
makePQ :: (Binary a, Ord a, Sizeable a) => PQ_Conf -> IO (PQ a)
makePQ _ = PQ `fmap` newIORef (Empty,0)

-- | Deletes the priority queue and all associated temporary files.
deletePQ :: PQ a -> IO ()
deletePQ (PQ _) = return ()

withPQ :: (Binary a, Ord a, Sizeable a) => PQ_Conf -> (PQ a -> IO b) -> IO b
withPQ conf = CE.bracket (makePQ conf) deletePQ

-- | Enqueues an element.
-- This operation may result in the creation of a file or in an enormous
-- merge of already created files.
enqueuePQ :: (Binary a, Ord a, Sizeable a) => a -> PQ a -> IO ()
enqueuePQ a (PQ pq) = do (p,s) <- readIORef pq
                         let !p' = insert a p
                             !s' = 1 + s
                         writeIORef pq (p',s')

-- | Removes the minimum element from the queue.
-- If the queue is already empty, nothing happens.  As a result, it is
-- possible that one or more file become empty and are deleted.
dequeuePQ :: (Binary a, Ord a, Sizeable a ) => PQ a -> IO ()
dequeuePQ (PQ pq) = do (p,s) <- readIORef pq
                       let !p' = dropMin p
                           !s' = max 0 (s - 1)
                       writeIORef pq (p',s')


-- | Returns the minimum element from the queue.
-- If the queue is empty, Nothing is returned.  Else the minimum element
-- currently in the queue.
peekMinPQ :: (Binary a, Ord a, Sizeable a) => PQ a -> IO (Maybe a)
peekMinPQ (PQ pq) = (getMin . fst) `fmap` readIORef pq

getMinPQ :: (Binary a, Ord a, Sizeable a) => PQ a -> IO (Maybe a)
getMinPQ (PQ pq) = do r <- (getMin . fst) `fmap` readIORef pq
                      case r of Nothing -> return () ; Just _ -> dequeuePQ  (PQ pq)
                      return r

sizePQ :: (Binary a, Ord a, Sizeable a) => PQ a -> IO Int
sizePQ (PQ pq) = snd `fmap` readIORef pq


-- We need an in-memory heap anyway.  Here's a skew heap.
data SkewHeap a = Empty | Node a (SkewHeap a) (SkewHeap a)

singleton :: Ord a => a -> SkewHeap a
singleton x = Node x Empty Empty

union :: Ord a => SkewHeap a -> SkewHeap a -> SkewHeap a
Empty              `union` t2                 = t2
t1                 `union` Empty              = t1
t1@(Node x1 l1 r1) `union` t2@(Node x2 l2 r2)
   | x1 <= x2                                 = Node x1 (t2 `union` r1) l1
   | otherwise                                = Node x2 (t1 `union` r2) l2

insert :: Ord a => a -> SkewHeap a -> SkewHeap a
insert x heap = singleton x `union` heap

getMin :: Ord a => SkewHeap a -> Maybe a
getMin Empty        = Nothing
getMin (Node x _ _) = Just x

dropMin :: Ord a => SkewHeap a -> SkewHeap a
dropMin Empty        = error "dropMin on empty queue... are you sure?!"
dropMin (Node _ l r) = l `union` r

