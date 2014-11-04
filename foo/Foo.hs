{-# LANGUAGE TemplateHaskell #-}
module Foo where

import Avro
import Language.Haskell.TH

data Foo = Foo { foo_one :: Int, bar_one :: String }
         | Bar { foo_two :: Int, bar_two :: String }
--                                  ^ field name?
--          ^ record type name?
--    ^ union type

data Bar = One | Two | Three

$( do reify ''Foo >>= runIO . print . ppr
      return [] )

deriveAvro ''Bar
deriveAvro ''Foo
