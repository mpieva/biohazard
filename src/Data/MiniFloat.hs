{-# LANGUAGE TypeFamilies, FlexibleInstances, CPP #-}
{-# LANGUAGE MultiParamTypeClasses, TemplateHaskell #-}
module Data.MiniFloat ( Mini(..), float2mini, mini2float ) where

import Data.Bits
import Data.Ix
import Data.Word                    ( Word8 )
import Data.Vector.Unboxed.Deriving ( derivingUnbox )
import Prelude

#if __GLASGOW_HASKELL__ == 704
import Data.Vector.Generic          ( Vector(..) )
import Data.Vector.Generic.Mutable  ( MVector(..) )
#endif

data Mini = Mini { unMini :: Word8 } deriving ( Eq, Ord, Show, Ix, Bounded )

derivingUnbox "Mini" [t| Mini -> Word8 |] [| unMini |] [| Mini |]

-- | Conversion to 0.4.4 format minifloat:  This minifloat fits into a
-- byte.  It has no sign, four bits of precision, and the range is from
-- 0 to 63488, initially in steps of 1/8.  Nice to store quality scores
-- with reasonable precision and range.
float2mini :: RealFloat a => a -> Mini
float2mini f | f' <  0   = error "no negative minifloats"   -- negative zero is fine!
             | f  <  2   = Mini f'
             | e >= 17   = Mini 0xff
             | s  < 16   = error $ "oops: " ++ show (e,s)
             | s  < 32   = Mini $ (e-1) `shiftL` 4 .|. (s .&. 0xf)
             | s == 32   = Mini $ e `shiftL` 4
             | otherwise = error $ "oops: " ++ show (e,s)
  where
    f' = round (8*f)
    e  = fromIntegral $ exponent f
    s  = round $ 32 * significand f

-- | Conversion from 0.4.4 format minifloat, see 'float2mini'.
mini2float :: Fractional a => Mini -> a
mini2float (Mini w) |  e == 0   =       fromIntegral w / 8.0
                    | otherwise = 2^e * fromIntegral m / 16.0
  where
    m = (w .&. 0xF) .|. 0x10
    e = w `shiftR` 4


