{-# LANGUAGE OverloadedStrings, TemplateHaskell #-}
-- Scan file with GT likelihoods, fit something...
--
-- First iteration:  Mitochondrion only.   We don't need to fit
-- anything.

import Bio.Bam.Pileup
import Bio.Genocall.AvroFile
import Bio.Iteratee
import Bio.Util
import Control.Monad ( zipWithM_, forM_ )
import Data.Avro
import Data.ByteString ( ByteString )
import Data.HashMap.Strict ( toList )
-- import Data.Iteratee
import Data.List ( foldl' )
import Data.Text.Encoding ( decodeUtf8 )
import Data.Text ( Text, unpack )

import qualified Data.Vector.Unboxed as U

{-
main = do -- let ctr :: ByteString
          ctr <- enumPure1Chunk test_data >=> run $ joinI $ writeAvroContainer ctr_opts $ stream2stream
          print ctr
          rd <- enumPure1Chunk (ctr :: ByteString) >=> run $ joinI $
                (readAvroContainer {- :: Enumeratee ByteString [Int] IO [Int] -}) $ stream2list

          print rd >> print (rd == test_data)

  where
    ctr_opts = ContainerOpts 1 "test"

    test_data :: [GenoCallBlock]
    test_data = [ GenoCallBlock "chr1" 0
                    [ GenoCallSite (CallStats 99 50 3500 99000)
                                   [0,1,2]
                                   (CallStats 80 40 3000 47000)
                                   [ IndelVariant 1 (V_Nuc $ fromList $ read "ACGT") ]
                                   [1001,1002,1003] ]
                , GenoCallBlock "MT" 0 [] ]
-}

main :: IO ()
main = enumDefaultInputs >=> run >=> print $
       joinI $ readAvroContainer $ \meta -> do
            -- liftIO . forM_ (toList meta) $ \(k,v) ->
                -- putStrLn $ unpack k ++ ": " ++ unpack (decodeUtf8 v)
            foldStream plus (0::Double)
  where
    plus acc cb = foldl' (\a cs -> a + U.minimum (U.map mini2float $ snp_likelihoods cs)) acc (called_sites cb)

pblock :: GenoCallBlock -> IO ()
pblock b = zipWithM_ psite [start_position b ..] (called_sites b)
  where
    psite p s = print (reference_name b, p, snp_likelihoods s)
