{-# LANGUAGE ForeignFunctionInterface #-}
module Bio.Util.MMap where

import Bio.Prelude
import Data.ByteString.Internal ( fromForeignPtr )
import Foreign.C.Types

unsafeMMapFile :: FilePath -> IO Bytes
unsafeMMapFile fp =
    bracket (openFd fp ReadOnly Nothing defaultFileFlags) closeFd $ \fd -> do
        stat <- getFdStatus fd
        let size = fromIntegral (fileSize stat)
        if size <= 0
            then return mempty
            else do
                ptr <- c_mmap size (fromIntegral fd)
                if ptr == nullPtr
                    then error "unable to mmap file"
                    else do
                          fptr <- newForeignPtrEnv c_munmap (intPtrToPtr $ fromIntegral size) ptr
                          return $ fromForeignPtr fptr 0 (fromIntegral size)

foreign import ccall unsafe  "my_mmap"   c_mmap   :: CSize -> CInt -> IO (Ptr Word8)
foreign import ccall unsafe "&my_munmap" c_munmap :: FunPtr (Ptr () -> Ptr Word8 -> IO ())


