{-# LANGUAGE CPP, TypeOperators, TypeSynonymInstances, FlexibleInstances #-}
module Bio.Prelude (
    module Bio.Base,
    module BasePrelude,
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

#if !MIN_VERSION_base(4,7,0)
    isLeft,
    isRight,
#endif

    Hashable(..),
    Unpack(..),
    hPutStr,
    hPutStrLn,
    fdPut,
    fdPutLazy,
    withFd,
    stderr,
    stdout,
    stdin
                   ) where

#if MIN_VERSION_base(4,7,0)
import BasePrelude  hiding ( EOF )
#else
import BasePrelude  hiding ( EOF, block )
#endif

import Bio.Base
import Data.ByteString     ( ByteString )
import Data.Text           ( Text )
import Data.Hashable       ( Hashable(..) )
import Data.HashMap.Strict ( HashMap )
import Data.HashSet        ( HashSet )
import Data.IntMap         ( IntMap )
import Data.IntSet         ( IntSet )
import Foreign.C.Error     ( throwErrnoIf_ )
import Foreign.Ptr         ( castPtr )
import System.IO           ( hPutStr, hPutStrLn, stderr, stdout, stdin )
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

#if !MIN_VERSION_base(4,7,0)
isLeft, isRight :: Either a b -> Bool
isLeft = either (const False) (const True)
isRight = either (const True) (const False)
#endif

fdPut :: Fd -> Bytes -> IO ()
fdPut fd s =
    B.unsafeUseAsCStringLen s $ \(p,l) ->
    throwErrnoIf_ (/= fromIntegral l) "fdPut" $
    fdWriteBuf fd (castPtr p) (fromIntegral l)

fdPutLazy :: Fd -> LazyBytes -> IO ()
fdPutLazy fd = mapM_ (fdPut fd) . BL.toChunks

withFd :: FilePath -> OpenMode -> Maybe FileMode -> OpenFileFlags
       -> (Fd -> IO a) -> IO a
withFd fp om fm ff k = bracket (openFd fp om fm ff) closeFd k
