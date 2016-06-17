{-# LANGUAGE CPP, TypeOperators, TypeSynonymInstances, FlexibleInstances #-}
module Bio.Prelude (
    module Bio.Base,
    module BasePrelude,
    Bytes,
    HashMap,
    HashSet,
    IntMap,
    IntSet,
    Text,
    Pair(..),
#ifndef __HADDOCK__
#ifdef __GLASGOW_HASKELL__
    (:!:),
#endif
#endif
    Hashable(..),
    Unpack(..),
    hPutStr,
    hPutStrLn,
    stderr,
    stdout,
    stdin
                   ) where

import Bio.Base
import BasePrelude  hiding ( EOF )

import Data.ByteString     ( ByteString )
import Data.Text           ( Text )
import Data.Hashable       ( Hashable(..) )
import Data.HashMap.Strict ( HashMap )
import Data.HashSet        ( HashSet )
import Data.IntMap         ( IntMap )
import Data.IntSet         ( IntSet )
import System.IO           ( hPutStr, hPutStrLn, stderr, stdout, stdin )

import qualified Data.ByteString.Char8 as S
import qualified Data.Text as T

type Bytes = ByteString

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
