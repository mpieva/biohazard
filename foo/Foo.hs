{-# LANGUAGE TemplateHaskell, OverloadedStrings #-}
import Bio.Iteratee
import Data.Avro
import Language.Haskell.TH
import qualified Data.Text as T
import System.IO
import Data.ByteString ( hPut ) 

data Foo = Foo { foo_one :: Int, bar_one :: T.Text }
         | Bar { foo_two :: Int, bar_two :: Double }
--                                  ^ field name?
--          ^ record type name?
--    ^ union type

data FooFoo = FooFoo { left :: Foo, right :: Foo }

data Bar = One | Two | Three

deriveAvro ''Bar
deriveAvro ''Foo
deriveAvro ''FooFoo

test_data = [
    FooFoo (Foo 23 "foo") (Bar 42 3.14),
    FooFoo (Foo 17 "und vier") (Foo 3  "bar") ]

main = withFile "foo.av" WriteMode $ \h ->
       enumPure1Chunk test_data >=> run $
       joinI $ writeAvroContainer (ContainerOpts 2 "experiment 1") $
       mapChunksM_ (hPut h)
