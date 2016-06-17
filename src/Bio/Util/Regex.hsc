{-# LANGUAGE CPP, ForeignFunctionInterface #-}
-- | The absolute minimum necessary for regex matching using POSIX regexec.
module Bio.Util.Regex ( Regex, regComp, regMatch ) where

#include <sys/types.h>
#include <regex.h>

import Foreign.Ptr
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc
import Foreign.C.String
import Foreign.C.Types
import System.IO.Unsafe
import BasePrelude

newtype Regex = Regex (ForeignPtr Regex)

regComp :: String -> Regex
regComp re = unsafePerformIO $ do
    fp <- mallocForeignPtrBytes #{size regex_t}
    withForeignPtr fp $ \p -> do
        withCString re $ \pre -> do
            ec <- regcomp p pre (#{const REG_EXTENDED} + #{const REG_NOSUB})
            when (ec /= 0) $ do
                sz <- regerror ec p nullPtr 0
                allocaBytes (fromIntegral sz) $ \err -> do
                    _ <- regerror ec p err sz
                    peekCString err >>= error . (++) "regexec: "
    addForeignPtrFinalizer regfree fp
    return $ Regex fp

regMatch :: Regex -> String -> Bool
regMatch (Regex fp) str =
    unsafePerformIO $
        withForeignPtr fp $ \p ->
            withCString str $ \s ->
                (==) 0 <$> regexec p s 0 nullPtr 0


foreign import ccall unsafe            regcomp :: Ptr Regex -> CString -> CInt -> IO CInt
foreign import ccall unsafe            regexec :: Ptr Regex -> CString -> CSize -> Ptr () -> CInt -> IO CInt
foreign import ccall unsafe           regerror :: CInt -> Ptr Regex -> CString -> CSize -> IO CSize
foreign import ccall unsafe "&regfree" regfree :: FunPtr (Ptr Regex -> IO ())
