{-# LANGUAGE FlexibleContexts, BangPatterns #-}

-- |Monadic Iteratees:
-- incremental input parsers, processors, and transformers
--
-- Iteratees for parsing binary data.

module Bio.Iteratee.Bytes (
  -- * Types
  Endian (..)
  -- * Endian multi-byte iteratees
  ,endianRead2
  ,endianRead3
  ,endianRead3i
  ,endianRead4
  ,endianRead8

  -- * Iteratees treating Bytes as list of Word8
  ,headStreamBS
  ,tryHeadBS
  ,peekStreamBS
  ,takeStreamBS
  ,dropStreamBS
  ,dropWhileStreamBS

  -- * Iteratees treating Bytes as list of Char
  ,enumLinesBS
  ,enumWordsBS
)
where

import Bio.Iteratee.Base
import Bio.Iteratee.Iteratee
import Bio.Prelude

import qualified Data.ByteString              as B
import qualified Data.ByteString.Char8        as C
import qualified Data.ByteString.Unsafe       as B

-- ------------------------------------------------------------------------
-- Binary Random IO Iteratees

-- Iteratees to read unsigned integers written in Big- or Little-endian ways

-- | Indicate endian-ness.
data Endian = MSB -- ^ Most Significant Byte is first (big-endian)
  | LSB           -- ^ Least Significan Byte is first (little-endian)
  deriving (Eq, Ord, Show, Enum)

endianRead2 :: Endian -> Iteratee Bytes m Word16
endianRead2 e = endianReadN e 2 word16'
{-# INLINE endianRead2 #-}

endianRead3 :: Endian -> Iteratee Bytes m Word32
endianRead3 e = endianReadN e 3 (word32' . (0:))
{-# INLINE endianRead3 #-}

-- |Read 3 bytes in an endian manner.  If the first bit is set (negative),
-- set the entire first byte so the Int32 will be negative as
-- well.
endianRead3i :: Monad m => Endian -> Iteratee Bytes m Int32
endianRead3i e = do
  c1 <- headStreamBS
  c2 <- headStreamBS
  c3 <- headStreamBS
  case e of
    MSB -> return $ (((fromIntegral c1
                        `shiftL` 8) .|. fromIntegral c2)
                        `shiftL` 8) .|. fromIntegral c3
    LSB ->
     let m :: Int32
         m = shiftR (shiftL (fromIntegral c3) 24) 8
     in return $ (((fromIntegral c3
                        `shiftL` 8) .|. fromIntegral c2)
                        `shiftL` 8) .|. fromIntegral m
{-# INLINE endianRead3i #-}

endianRead4 :: Endian -> Iteratee Bytes m Word32
endianRead4 e = endianReadN e 4 word32'
{-# INLINE endianRead4 #-}

endianRead8 :: Endian -> Iteratee Bytes m Word64
endianRead8 e = endianReadN e 8 word64'
{-# INLINE endianRead8 #-}

-- This function does all the parsing work, depending upon provided arguments
endianReadN ::
  Endian
  -> Int
  -> ([Word8] -> b)
  -> Iteratee Bytes m b
endianReadN MSB n0 cnct = liftI (step n0 [])
 where
  step !n acc (Chunk c)
    | B.null c        = liftI (step n acc)
    | B.length c >= n = let (this,next) = B.splitAt n c
                            !result     = cnct $ acc ++ B.unpack this
                        in idone result (Chunk next)
    | otherwise        = liftI (step (n - B.length c) (acc ++ B.unpack c))
  step !n acc (EOF Nothing)  = icont (step n acc) (Just $ toException EofException)
  step !n acc (EOF (Just e)) = icont (step n acc) (Just e)
endianReadN LSB n0 cnct = liftI (step n0 [])
 where
  step !n acc (Chunk c)
    | B.null c        = liftI (step n acc)
    | B.length c >= n = let (this,next) = B.splitAt n c
                            !result = cnct $ B.unpack (B.reverse this) ++ acc
                        in idone result (Chunk next)
    | otherwise        = liftI (step (n - B.length c)
                                     (B.unpack (B.reverse c) ++ acc))
  step !n acc (EOF Nothing)  = icont (step n acc)
                                    (Just $ toException EofException)
  step !n acc (EOF (Just e)) = icont (step n acc) (Just e)
{-# INLINE endianReadN #-}


word16' :: [Word8] -> Word16
word16' [c1,c2] = word16 c1 c2
word16' _ = error "iteratee: internal error in word16'"

word16 :: Word8 -> Word8 -> Word16
word16 c1 c2 = (fromIntegral c1 `shiftL`  8) .|.  fromIntegral c2
{-# INLINE word16 #-}

word32' :: [Word8] -> Word32
word32' [c1,c2,c3,c4] = word32 c1 c2 c3 c4
word32' _ = error "iteratee: internal error in word32'"

word32 :: Word8 -> Word8 -> Word8 -> Word8 -> Word32
word32 c1 c2 c3 c4 =
  (fromIntegral c1 `shiftL` 24) .|.
  (fromIntegral c2 `shiftL` 16) .|.
  (fromIntegral c3 `shiftL`  8) .|.
   fromIntegral c4
{-# INLINE word32 #-}

word64' :: [Word8] -> Word64
word64' [c1,c2,c3,c4,c5,c6,c7,c8] = word64 c1 c2 c3 c4 c5 c6 c7 c8
word64' _ = error "iteratee: internal error in word64'"
{-# INLINE word64' #-}

word64
  :: Word8 -> Word8 -> Word8 -> Word8
  -> Word8 -> Word8 -> Word8 -> Word8
  -> Word64
word64 c1 c2 c3 c4 c5 c6 c7 c8 =
  (fromIntegral c1 `shiftL` 56) .|.
  (fromIntegral c2 `shiftL` 48) .|.
  (fromIntegral c3 `shiftL` 40) .|.
  (fromIntegral c4 `shiftL` 32) .|.
  (fromIntegral c5 `shiftL` 24) .|.
  (fromIntegral c6 `shiftL` 16) .|.
  (fromIntegral c7 `shiftL`  8) .|.
   fromIntegral c8
{-# INLINE word64 #-}

headStreamBS :: Iteratee Bytes m Word8
headStreamBS = liftI step
  where
  step (Chunk c)
    | B.null c = icont step Nothing
    | otherwise  = idone (B.unsafeHead c) (Chunk (B.unsafeTail c))
  step stream          = icont step (Just (setEOF stream))
{-# INLINE headStreamBS #-}

peekStreamBS :: Iteratee Bytes m (Maybe Word8)
peekStreamBS = liftI step
  where
    step s@(Chunk vec)
      | B.null vec = liftI step
      | otherwise   = idone (Just $ B.unsafeHead vec) s
    step stream     = idone Nothing stream
{-# INLINE peekStreamBS #-}

tryHeadBS :: Iteratee Bytes m (Maybe Word8)
tryHeadBS = liftI step
  where
  step (Chunk vec)
    | B.null vec = liftI step
    | otherwise  = idone (Just (B.unsafeHead vec)) (Chunk (B.unsafeTail vec))
  step stream          = idone Nothing stream
{-# INLINE tryHeadBS #-}

dropStreamBS :: Int -> Iteratee Bytes m ()
dropStreamBS 0  = idone () (Chunk emptyP)
dropStreamBS n' = liftI (step n')
  where
    step n (Chunk str)
      | B.length str < n = liftI (step (n - B.length str))
      | otherwise        = idone () (Chunk (B.drop n str))
    step _ stream        = idone () stream
{-# INLINE dropStreamBS #-}

dropWhileStreamBS :: (Word8 -> Bool) -> Iteratee Bytes m ()
dropWhileStreamBS p = liftI step
  where
    step (Chunk str)
      | B.null rest  = liftI step
      | otherwise    = idone () (Chunk rest)
      where
        rest = B.dropWhile p str
    step stream      = idone () stream
{-# INLINE dropWhileStreamBS #-}

takeStreamBS ::
  Monad m
  => Int   -- ^ number of elements to consume
  -> Enumeratee Bytes Bytes m a
takeStreamBS n' iter
 | n' <= 0   = return iter
 | otherwise = Iteratee $ \od oc -> runIter iter (on_done od oc) (on_cont od oc)
  where
    on_done od oc x _ = runIter (dropStreamBS n' >> return (return x)) od oc
    on_cont od oc k Nothing = if n' == 0 then od (liftI k) (Chunk mempty)
                                 else runIter (liftI (step n' k)) od oc
    on_cont od oc _ (Just e) = runIter (dropStreamBS n' >> throwErr e) od oc
    step n k (Chunk str)
      | B.null str        = liftI (step n k)
      | B.length str <= n = takeStreamBS (n - B.length str) $ k (Chunk str)
      | otherwise          = idone (k (Chunk s1)) (Chunk s2)
      where (s1, s2) = B.splitAt n str
    step _n k stream       = idone (liftI k) stream
{-# INLINE takeStreamBS #-}

-- Like enumWords, but operates on ByteStrings.
-- This is provided as a higher-performance alternative to enumWords, and
-- is equivalent to treating the stream as a Data.ByteString.Char8.ByteString.
enumWordsBS :: Monad m => Enumeratee Bytes [Bytes] m a
enumWordsBS iter = convStream getter iter
  where
    getter = liftI step
    lChar = isSpace . C.last
    step (Chunk xs)
      | C.null xs  = getter
      | lChar xs   = idone (C.words xs) (Chunk C.empty)
      | otherwise  = icont (step' xs) Nothing
    step str       = idone mempty str
    step' xs (Chunk ys)
      | C.null ys  = icont (step' xs) Nothing
      | lChar ys   = idone (C.words . C.append xs $ ys) mempty
      | otherwise  = let w' = C.words . C.append xs $ ys
                         ws = init w'
                         ck = last w'
                     in idone ws (Chunk ck)
    step' xs str   = idone (C.words xs) str
{-# INLINE enumWordsBS #-}

-- Like enumLines, but operates on ByteStrings.
-- This is provided as a higher-performance alternative to enumLines, and
-- is equivalent to treating the stream as a Data.ByteString.Char8.ByteString.
enumLinesBS :: Monad m => Enumeratee Bytes [Bytes] m a
enumLinesBS = convStream getter
  where
    getter = icont step Nothing
    lChar = (== '\n') . C.last
    step (Chunk xs)
      | C.null xs  = getter
      | lChar xs   = idone (C.lines xs) (Chunk C.empty)
      | otherwise  = icont (step' xs) Nothing
    step str       = idone mempty str
    step' xs (Chunk ys)
      | C.null ys  = icont (step' xs) Nothing
      | lChar ys   = idone (C.lines . C.append xs $ ys) mempty
      | otherwise  = let w' = C.lines $ C.append xs ys
                         ws = init w'
                         ck = last w'
                     in idone ws (Chunk ck)
    step' xs str   = idone (C.lines xs) str
{-# INLINE enumLinesBS #-}
