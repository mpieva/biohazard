{-# LANGUAGE TemplateHaskell #-}
module Foo where

import Avro
import Language.Haskell.TH
import qualified Data.Text as T

data Foo = Foo { foo_one :: Int, bar_one :: T.Text }
         | Bar { foo_two :: Int, bar_two :: Double }
--                                  ^ field name?
--          ^ record type name?
--    ^ union type

data FooFoo = FooFoo { left :: Foo, right :: Foo }

data Bar = One | Two | Three

$( do reify ''Foo >>= runIO . print . ppr
      return [] )

deriveAvro ''Bar
deriveAvro ''Foo
deriveAvro ''FooFoo
