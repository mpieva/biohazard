{-# LANGUAGE CPP, TypeOperators #-}
module Bio.Prelude (
    module Bio.Base,
    module BasePrelude,
    module System.IO,
    module System.Posix.Files,
    module System.Posix.IO,
    module System.Posix.Types,

    Bytes, LazyBytes,
    HashMap,
    HashSet,
    IntMap,
    IntSet,
    Text, LazyText,
    Pair(..),
#ifndef __HADDOCK__
#ifdef __GLASGOW_HASKELL__
    (:!:),
#endif
#endif

#if !MIN_VERSION_base(4,8,0)
    first,
    second,
#endif

    decodeBytes,
    encodeBytes,

    Hashable(..),
    Unpack(..),
    fdGet,
    fdPut,
    fdPutLazy,
    withFd
                   ) where

import BasePrelude
#if MIN_VERSION_base(4,9,0)
                    hiding ( EOF, log1p, log1pexp, log1mexp, expm1 )
#else
                    hiding ( EOF )
#endif

#if !MIN_VERSION_base(4,8,0)
-- Not as nice as Data.Bifunctor, but still useful.
import Control.Arrow       ( first, second )
#endif

import Bio.Base
import Data.ByteString     ( ByteString )
import Data.ByteString.Internal ( createAndTrim )
import Data.Text           ( Text )
import Data.Hashable       ( Hashable(..) )
import Data.HashMap.Strict ( HashMap )
import Data.HashSet        ( HashSet )
import Data.IntMap         ( IntMap )
import Data.IntSet         ( IntSet )
import Data.Text.Encoding  ( encodeUtf8, decodeUtf8With )
import Foreign.C.Error     ( throwErrnoIf_ )
import System.IO           ( hPrint, hPutStr, hPutStrLn, stderr, stdout, stdin )
import System.Posix.Files
import System.Posix.IO
import System.Posix.Types

import qualified Data.ByteString.Unsafe as B
import qualified Data.ByteString.Lazy   as BL
import qualified Data.ByteString.Char8  as S
import qualified Data.Text              as T
import qualified Data.Text.Lazy         as TL

type Bytes     =    ByteString
type LazyBytes = BL.ByteString
type LazyText  = TL.Text

infixl 2 :!:

-- | A strict pair.
data Pair a b = !a :!: !b deriving(Eq, Ord, Show, Read, Bounded, Ix)

#ifndef __HADDOCK__
#ifdef __GLASGOW_HASKELL__
-- This gives a nicer syntax for the type but only works in GHC for now.
type (:!:) = Pair
#endif
#endif

-- | Class of things that can be unpacked into 'String's.  Kind of the
-- opposite of 'IsString'.
class Unpack s where unpack :: s -> String

instance Unpack ByteString where unpack = S.unpack
instance Unpack Text       where unpack = T.unpack
instance Unpack String     where unpack = id

-- | @fdGet bs fd@ reads up to @bs@ 'Bytes' from file descriptor @Fd@.
-- Returns an empty 'Bytes' at end of file.
fdGet :: Int -> Fd -> IO Bytes
fdGet bs fd =
    createAndTrim bs $ \p ->
        fromIntegral <$> fdReadBuf fd (castPtr p) (fromIntegral bs)

fdPut :: Fd -> Bytes -> IO ()
fdPut fd s = B.unsafeUseAsCStringLen s $ \(p,l) ->
             throwErrnoIf_ (/= fromIntegral l) "fdPut" $
             fdWriteBuf fd (castPtr p) (fromIntegral l)

fdPutLazy :: Fd -> LazyBytes -> IO ()
fdPutLazy fd = mapM_ (fdPut fd) . BL.toChunks

withFd :: FilePath -> OpenMode -> Maybe FileMode -> OpenFileFlags
       -> (Fd -> IO a) -> IO a
withFd fp om fm ff = bracket (openFd fp om fm ff) closeFd

-- | Converts 'Bytes' into 'Text'.  This uses UTF8, but if there is an
-- error, it pretends it was Latin1.  Evil as this is, it tends to Just
-- Work on files where nobody ever wasted a thought on encodings.
decodeBytes :: Bytes -> Text
decodeBytes = decodeUtf8With (const $ fmap w2c)

-- | Converts 'Text' into 'Bytes'.  This uses UTF8.
encodeBytes :: Text -> Bytes
encodeBytes = encodeUtf8

