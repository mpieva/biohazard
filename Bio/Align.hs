{-# LANGUAGE ForeignFunctionInterface #-}
module Bio.Align (
    Mode(..),
    myersAlign,
    showAligned
                 ) where

import Control.Applicative
import Foreign.C.String
import Foreign.C.Types
import Foreign.Marshal.Alloc
import System.IO.Unsafe ( unsafePerformIO )

import qualified Data.ByteString.Char8      as S
import qualified Data.ByteString.Unsafe     as S
import qualified Data.ByteString.Lazy.Char8 as L

foreign import ccall unsafe "myers_align.h myers_diff" myers_diff ::
        CString -> CInt ->              -- sequence A and length A
        CInt ->                         -- mode (an enum)
        CString -> CInt ->              -- sequence B and length B
        CInt ->                         -- max distance
        CString ->                      -- backtracing space A
        CString ->                      -- backtracing space B
        IO CInt                         -- returns distance

data Mode = Globally | IsPrefix | HasPrefix deriving Enum

myersAlign :: S.ByteString -> S.ByteString -> Mode -> Int -> (Int, S.ByteString, S.ByteString)
myersAlign seqA seqB mode maxd =
    unsafePerformIO                                 $
    S.unsafeUseAsCStringLen seqA                    $ \(seq_a, len_a) ->
    S.unsafeUseAsCStringLen seqB                    $ \(seq_b, len_b) ->

    -- size of output buffers derives from this:
    -- char *out_a = bt_a + len_a + maxd +2 ;
    -- char *out_b = bt_b + len_b + maxd +2 ;
    allocaBytes (len_a + maxd + 2)                  $ \bt_a ->
    allocaBytes (len_b + maxd + 2)                  $ \bt_b ->

    myers_diff seq_a (fromIntegral len_a)
               (fromIntegral $ fromEnum mode)
               seq_b (fromIntegral len_b)
               (fromIntegral maxd) bt_a bt_b      >>= \dist ->
    (,,) (fromIntegral dist) <$>
         S.packCString bt_a  <*>
         S.packCString bt_b

showAligned :: Int -> [S.ByteString] -> L.ByteString
showAligned w ss | all S.null ss = L.empty
                 | otherwise = L.concat [ L.unlines $ map (L.fromChunks . (:[])) lefts
                                        , L.pack $ agreement ++ "\n\n"
                                        , showAligned w rights ]
  where
    (lefts, rights) = unzip $ map (S.splitAt w) ss
    agreement = map star $ S.transpose lefts
    star str = if S.null str || S.all (== S.head str) str then ' ' else '*'

