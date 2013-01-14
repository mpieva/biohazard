module Bio.PriorityQueue (
        Sizable(..),
        PQ,
        makePQ,
        deletePQ,
        enqueuePQ,
        dequeuePQ,
        minimumPQ
        ) where

-- | A Priority Queue that can fall back to external storage.  
-- 
-- Whatever is to be stored in this queue needs to be in Binary, because
-- it may need to be moved to external storage on demand.  We also need
-- a way to estimate the memory consumption of an enqueued object.  When
-- constructing the queue, two limits are set, (i) the maximum amount of
-- RAM to consume, and (ii) the maximum number of files to keep open.
-- Note that open files use memory for buffering, too.
--
-- Enqueued objects are kept in an in memory heap until the memory
-- consumption becomes too high.  At that point, a file is created, the
-- whole heap is sorted and dumped, then the file is reopened.  Open
-- files are again kept in a heap so that dequeueing objects amounts to
-- performing a merge sort on multiple files.  If too many files are
-- open, we do exactly that:  merge-sort all files and the in-memory
-- heap into a single new file. 
--
-- XXX  Truth be told, this queue isn't backed externally, and ignores
--      all limits.  It *is* a Priority Queue, though!
--
-- XXX  May want to add callbacks for significant events (new file,
--      massive merge, deletion of file?)

import Data.Binary
import Data.IORef

newtype PQ a = PQ (IORef (SkewHeap a))

newtype MaxMb = MaxMb Int
newtype MaxFd = MaxFd Int

class Sizable a where usedBytes :: a -> Int

-- | Creates a priority queue.  Note that the priority queue creates
-- files, which will only be cleaned up if deletePQ is called.
makePQ :: (Binary a, Ord a, Sizable a)
       => MaxMb             -- ^ memory limit
       -> MaxFd             -- ^ max number of open files
       -> FilePath          -- ^ path to temporary files (a directory will be created)
       -> IO (PQ a)         -- ^ makes a priority queue
makePQ (MaxMb _) (MaxFd _) _tmp = PQ `fmap` newIORef Empty

-- | Deletes the priority queue and all associated temporary files.
--
deletePQ :: PQ a -> IO ()
deletePQ (PQ _) = return ()

-- | Enqueues an element.
-- This operation may result in the creation of a file or in an enormous
-- merge of already created files.
enqueuePQ :: (Binary a, Ord a, Sizable a) => a -> PQ a -> IO ()
enqueuePQ a (PQ pq) = do p <- readIORef pq
                         writeIORef pq $! insert a p

-- | Removes the minimum element from the queue.
-- If the queue is already empty, nothing happens.  As a result, it is
-- possible that one or more file become empty and are deleted.
dequeuePQ :: (Binary a, Ord a, Sizable a) => PQ a -> IO ()
dequeuePQ (PQ pq) = do p <- readIORef pq
                       writeIORef pq $! dropMin p


-- | Returns the minimum element from the queue.  
-- If the queue is empty, Nothing is returned.  Else the minimum element
-- currently in the queue.
minimumPQ :: (Binary a, Ord a, Sizable a) => PQ a -> IO (Maybe a)
minimumPQ (PQ pq) = getMin `fmap` readIORef pq


-- We need an in-memory heap anyway.  Here's a skew heap.
data SkewHeap a = Empty
                | Node a (SkewHeap a) (SkewHeap a)
 
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
dropMin Empty        = Empty
dropMin (Node _ l r) = l `union` r
