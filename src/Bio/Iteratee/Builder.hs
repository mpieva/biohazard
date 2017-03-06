{-# LANGUAGE ForeignFunctionInterface #-}
-- | Buffer builder to assemble Bgzf blocks.  The plan is to serialize
-- stuff (BAM and BCF) into a buffer, then Bgzf chunks from the buffer.
-- We use a large buffer, and we always make sure there is plenty of
-- space in it (to avoid redundant checks).  Whenever a block is ready
-- to be compressed, we stick it into a MVar.  When we run out of space,
-- we simply use a new buffer.  Multiple threads grab pieces from the
-- MVar, compress them, pass them downstream through another MVar.  A
-- final thread restores the order and writes the blocks.

module Bio.Iteratee.Builder (
    BB(..),
    newBuffer,
    fillBuffer,
    expandBuffer,
    encodeBgzf,
    BgzfTokens(..),
    BclArgs(..),
    BclSpecialType(..),
    int_loop,
    loop_bcl_special
                            ) where

import Bio.Iteratee
import Bio.Iteratee.Bgzf                   ( compressChunk, maxBlockSize, bgzfEofMarker )
import Bio.Prelude
import Foreign.ForeignPtr                  ( ForeignPtr, withForeignPtr, mallocForeignPtrBytes )
import Foreign.Marshal.Utils               ( copyBytes )
import Foreign.Ptr                         ( Ptr, plusPtr )
import Foreign.Storable                    ( pokeByteOff )

import qualified Data.ByteString            as B
import qualified Data.ByteString.Unsafe     as B
import qualified Data.Vector.Storable       as VS

-- | We manage a large buffer (multiple megabytes), of which we fill an
-- initial portion.  We remeber the size, the used part, and two marks
-- where we later fill in sizes for the length prefixed BAM or BCF
-- records.  We move the buffer down when we yield a piece downstream,
-- and when we run out of space, we simply move to a new buffer.
-- Garbage collection should take care of the rest.  Unused 'mark' must
-- be set to (maxBound::Int) so it doesn't interfere with flushing.

data BB = BB { buffer :: {-# UNPACK #-} !(ForeignPtr Word8)
             , size   :: {-# UNPACK #-} !Int            -- total size of buffer
             , off    :: {-# UNPACK #-} !Int            -- offset of active portion
             , used   :: {-# UNPACK #-} !Int            -- used portion (inactive & active)
             , mark   :: {-# UNPACK #-} !Int            -- offset of mark
             , mark2  :: {-# UNPACK #-} !Int }          -- offset of mark2

instance Show BB where
    show bb = show (size bb, off bb, used bb, mark bb, mark2 bb)

-- | Things we are able to encode.  Taking inspiration from
-- binary-serialise-cbor, we define these as a lazy list-like thing and
-- consume it in a interpreter.

data BgzfTokens = TkWord32   {-# UNPACK #-} !Word32       BgzfTokens -- a 4-byte int
                | TkWord16   {-# UNPACK #-} !Word16       BgzfTokens -- a 2-byte int
                | TkWord8    {-# UNPACK #-} !Word8        BgzfTokens -- a byte
                | TkFloat    {-# UNPACK #-} !Float        BgzfTokens -- a float
                | TkDouble   {-# UNPACK #-} !Double       BgzfTokens -- a double
                | TkString   {-# UNPACK #-} !B.ByteString BgzfTokens -- a raw string
                | TkDecimal  {-# UNPACK #-} !Int          BgzfTokens -- roughly ':%d'
                | TkLnString {-# UNPACK #-} !B.ByteString BgzfTokens -- a length-prefixed string
                -- lotsa stuff might be missing here

                | TkSetMark                               BgzfTokens -- sets the first mark
                | TkEndRecord                             BgzfTokens -- completes a BAM record
                | TkEndRecordPart1                        BgzfTokens -- completes part 1 of a BCF record
                | TkEndRecordPart2                        BgzfTokens -- completes part 2 of a BCF record
                | TkEnd                                              -- nothing more, for now

                -- specialties
                | TkBclSpecial !BclArgs                   BgzfTokens
                | TkLowLevel {-# UNPACK #-} !Int (BB -> IO BB) BgzfTokens

data BclSpecialType = BclNucsBin | BclNucsAsc | BclNucsAscRev | BclQualsBin | BclQualsAsc | BclQualsAscRev

data BclArgs = BclArgs BclSpecialType
                       {-# UNPACK #-} !(VS.Vector Word8)  -- bcl matrix
                       {-# UNPACK #-} !Int                -- stride
                       {-# UNPACK #-} !Int                -- first cycle
                       {-# UNPACK #-} !Int                -- last cycle
                       {-# UNPACK #-} !Int                -- cluster index

-- | Creates a buffer.
newBuffer :: Int -> IO BB
newBuffer sz = mallocForeignPtrBytes sz >>= \ar -> return $ BB ar sz 0 0 maxBound maxBound

-- | Creates a new buffer, copying the content from an old one, with
-- higher capacity.
expandBuffer :: Int -> BB -> IO BB
expandBuffer minsz b = do
    let sz' = max (2 * (size b - used b)) minsz
    arr1 <- mallocForeignPtrBytes sz'
    withForeignPtr arr1 $ \d ->
        withForeignPtr (buffer b) $ \s ->
             copyBytes d (plusPtr s (off b)) (used b - off b)
    return $ BB { buffer = arr1
                , size   = sz'
                , off    = 0
                , used   = used b - off b
                , mark   = if mark  b == maxBound then maxBound else mark  b - off b
                , mark2  = if mark2 b == maxBound then maxBound else mark2 b - off b }

compressChunk' :: Int -> ForeignPtr Word8 -> Int -> Int -> IO B.ByteString
compressChunk' lv fptr off len =
    withForeignPtr fptr $ \ptr ->
    compressChunk lv (plusPtr ptr off) (fromIntegral len)

instance Nullable (Endo BgzfTokens) where
    nullC f = case appEndo f TkEnd of TkEnd -> True ; _ -> False

-- | Expand a chain of tokens into a buffer, sending finished pieces
-- downstream as soon as possible.
encodeBgzf :: MonadIO m => Int -> Enumeratee (Endo BgzfTokens) B.ByteString m b
encodeBgzf lv = (\out -> newBuffer (1024*1024) `ioBind` \bb -> eneeCheckIfDone (liftI . go bb) out)
                ><> parRunIO (2*numCapabilities)
  where
    go bb0 k (EOF  mx) = final_flush bb0 mx k
    go bb0 k (Chunk f)
        -- initially, we make sure we have reasonable space.  this may not be enough.
        | size bb0 - used bb0 < 1024 = expandBuffer (1024*1024) bb0 `ioBind` \bb' -> go' bb' k (appEndo f TkEnd)
        | otherwise                  =                                               go' bb0 k (appEndo f TkEnd)

    -- we arrive here because we ran out of buffer space, so we always expand it.
    go1 bb0 k tk = expandBuffer (1024*1024) bb0 `ioBind` \bb' -> go' bb' k tk

    go' bb0 k tk = fillBuffer bb0 tk `ioBind` \(bb',tk') -> flush_blocks tk' bb' k


    -- We can flush anything that is between 'off' and the lower of 'mark'
    -- and 'used'.  When done, we bump 'off'.
    flush_blocks tk bb k
        | min (mark bb) (used bb) - off bb < maxBlockSize =
            case tk of TkEnd -> liftI $ go bb k
                       _     -> go1 bb k tk

        | otherwise = do
            eneeCheckIfDone (flush_blocks tk bb { off = off bb + maxBlockSize }) $
                k $ Chunk [compressChunk' lv (buffer bb) (off bb) maxBlockSize]

    final_flush bb mx k
        | used bb > off bb =
            idone (k $ Chunk [ compressChunk' lv (buffer bb) (off bb) (used bb - off bb)
                             , return bgzfEofMarker ]) (EOF mx)
        | otherwise =
            idone (k $ Chunk [ return bgzfEofMarker ]) (EOF mx)


fillBuffer :: BB -> BgzfTokens -> IO (BB, BgzfTokens)
fillBuffer bb0 tk = withForeignPtr (buffer bb0) (\p -> go_slowish p bb0 tk)
  where
    go_slowish p bb tk1 = go_fast p bb (used bb) tk1

    go_fast p bb use tk1 = case tk1 of
        -- no space?  not our job.
        _ | size bb - use < 1024 -> return (bb { used = use },tk1)

        -- the actual end.
        TkEnd                    -> return (bb { used = use },tk1)

        -- I'm cheating.  This stuff works only if the platform allows
        -- unaligned accesses, is little-endian and uses IEEE floats.
        -- It's true on i386 and ix86_64.
        TkWord32   x tk' -> do pokeByteOff p use x
                               go_fast p bb (use + 4) tk'

        TkWord16   x tk' -> do pokeByteOff p use x
                               go_fast p bb (use + 2) tk'

        TkWord8    x tk' -> do pokeByteOff p use x
                               go_fast p bb (use + 1) tk'

        TkFloat    x tk' -> do pokeByteOff p use x
                               go_fast p bb (use + 4) tk'

        TkDouble   x tk' -> do pokeByteOff p use x
                               go_fast p bb (use + 8) tk'

        TkString   s tk'
            -- Too big, can't handle.  We will get progressively bigger
            -- buffers and eventually handle it; for very large strings,
            -- it works, but isn't ideal.  XXX
            | B.length s > size bb - use -> return (bb { used = use },tk')

            | otherwise  -> do let ln = B.length s
                               B.unsafeUseAsCString s $ \q ->
                                    copyBytes (p `plusPtr` use) q ln
                               go_fast p bb (use + ln) tk'

        TkDecimal  x tk' -> do ln <- int_loop (p `plusPtr` use) x
                               go_fast p bb (use + ln) tk'

        TkLnString s tk'
            -- Too big, can't handle.  We will get progressively bigger
            -- buffers and eventually handle it; for very large strings,
            -- it works, but isn't ideal.  XXX
            | B.length s > size bb - use - 4 -> return (bb { used = use },tk')

            | otherwise  -> do let ln = B.length s
                               pokeByteOff p use (fromIntegral ln :: Word32)
                               B.unsafeUseAsCString s $ \q ->
                                    copyBytes (p `plusPtr` (use + 4)) q ln
                               go_fast p bb (use + ln + 4) tk'

        TkSetMark        tk' ->    go_slowish p bb { used = use + 4, mark = use } tk'

        TkEndRecord      tk' -> do let !l = use - mark bb - 4
                                   pokeByteOff p (mark bb) (fromIntegral l :: Word32)
                                   go_slowish p bb { used = use, mark = maxBound } tk'

        TkEndRecordPart1 tk' -> do let !l = use - mark bb - 4
                                   pokeByteOff p (mark bb - 4) (fromIntegral l :: Word32)
                                   go_slowish p bb { used = use, mark2 = use } tk'

        TkEndRecordPart2 tk' -> do let !l = use - mark2 bb
                                   pokeByteOff p (mark bb) (fromIntegral l :: Word32)
                                   go_slowish p bb { used = use, mark = maxBound } tk'


        TkBclSpecial special_args tk' -> do
            l <- loop_bcl_special (p `plusPtr` use) special_args
            go_fast p bb (use + l) tk'

        TkLowLevel minsize proc tk'
            | size bb - use < minsize -> return (bb { used = use },tk1)
            | otherwise               -> flip (,) tk' <$> proc bb


loop_bcl_special :: Ptr Word8 -> BclArgs -> IO Int
loop_bcl_special p (BclArgs tp vec stride u v i) =

    VS.unsafeWith vec $ \q -> case tp of
        BclNucsBin -> do
            nuc_loop p stride (plusPtr q i) u v
            return $ (v - u + 2) `div` 2

        BclNucsAsc -> do
            nuc_loop_asc p stride (plusPtr q i) u v
            return $ v - u + 1

        BclNucsAscRev -> do
            nuc_loop_asc_rev p stride (plusPtr q i) u v
            return $ v - u + 1

        BclQualsBin -> do
            qual_loop p stride (plusPtr q i) u v
            return $ v - u + 1

        BclQualsAsc -> do
            qual_loop_asc p stride (plusPtr q i) u v
            return $ v - u + 1

        BclQualsAscRev -> do
            qual_loop_asc_rev p stride (plusPtr q i) u v
            return $ v - u + 1

foreign import ccall unsafe "nuc_loop"
    nuc_loop :: Ptr Word8 -> Int -> Ptr Word8 -> Int -> Int -> IO ()

foreign import ccall unsafe "nuc_loop_asc"
    nuc_loop_asc :: Ptr Word8 -> Int -> Ptr Word8 -> Int -> Int -> IO ()

foreign import ccall unsafe "nuc_loop_asc_rev"
    nuc_loop_asc_rev :: Ptr Word8 -> Int -> Ptr Word8 -> Int -> Int -> IO ()

foreign import ccall unsafe "qual_loop"
    qual_loop :: Ptr Word8 -> Int -> Ptr Word8 -> Int -> Int -> IO ()

foreign import ccall unsafe "qual_loop_asc"
    qual_loop_asc :: Ptr Word8 -> Int -> Ptr Word8 -> Int -> Int -> IO ()

foreign import ccall unsafe "qual_loop_asc_rev"
    qual_loop_asc_rev :: Ptr Word8 -> Int -> Ptr Word8 -> Int -> Int -> IO ()

foreign import ccall unsafe "int_loop"
    int_loop :: Ptr Word8 -> Int -> IO Int

