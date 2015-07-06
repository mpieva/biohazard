{-# LANGUAGE TemplateHaskell, OverloadedStrings, BangPatterns     #-}
{-# LANGUAGE MultiParamTypeClasses, TypeFamilies, RecordWildCards #-}
{-# LANGUAGE ForeignFunctionInterface, GeneralizedNewtypeDeriving #-}

-- Two-stage demultiplexing.
--
-- We assume we know the list of i7 and i5 index oligos.  We seek to
-- decompose a set of reads into a mix of pairs of these by the Maximum
-- Likelihood method.  Once that's done, an empirical Bayesian Maximum
-- Posterior call is done.  All kinds of errors can be rolled into one
-- quality score.
--
-- TODO
--
--  - Input layer to gather index sequences.  (Got a lame version.)
--  - First pass to gather data.  Any index read shall be represented
--    in a single Word64.  (Did that.  It's slow.)
--  - Multiple passes of the EM algorithm.  We start with a mix of
--    mostly known pairs and a little bit of all other pairs.
--  - Final calling pass from BAM to BAM.
--  - Auxillary statistics:  composition of the mix, false assignment
--    matrix or rates per read group, maximum achievable false
--    assignment rates.

import Bio.Bam
import Control.Arrow ( (&&&) )
import Control.Monad ( forM_, when )
import Control.Monad.Fix
import Control.Monad.ST
import Control.Monad.IO.Class
import Data.Aeson
import Data.Bits
import Data.Char ( chr )
import Data.List ( foldl' )
import Data.Vector.Unboxed.Deriving
import Data.Word ( Word64 )
import Foreign.Ptr
import Foreign.Storable
import System.Directory ( getHomeDirectory )
import System.IO
import System.IO.Unsafe
import System.Random ( randomRIO )
import System.Time ( getClockTime )

import qualified Data.ByteString as B
import qualified Data.HashMap.Strict as HM
import qualified Data.Text as T
import qualified Data.Text.Encoding as T
import qualified Data.Text.Format as T
import qualified Data.Vector as V
import qualified Data.Vector.Algorithms.Intro as V
import qualified Data.Vector.Unboxed as U
-- import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import qualified Data.Vector.Generic            as VG
import qualified Data.Vector.Generic.Mutable    as VGM


-- | An index sequence must have at most eight bases.  We represent a
-- base and its quality score in a single byte:  the top three bits are
-- the base ("ACGTN" = [0,1,3,2,7]), the lower five bits are the quality,
-- clamped to 31.

newtype Index = Index Word64 deriving Storable

instance Show Index where
    show (Index x) = [ "ACTGNNNN" !! fromIntegral b | i <- [56,48..0], let b = (x `shiftR` (i+5)) .&. 0x7 ]
            ++ 'q' : [ chr (fromIntegral q+33)      | i <- [56,48..0], let q = (x `shiftR` i) .&. 0x1F ]

derivingUnbox "Index" [t| Index -> Word64 |] [| \ (Index i) -> i |] [| Index |]

a `roundup` b = ((a + b - 1) `div` b) * b

fromS :: B.ByteString -> Index
fromS sq = fromSQ sq (B.replicate (B.length sq) 33)

fromSQ :: B.ByteString -> B.ByteString -> Index
fromSQ sq qs = Index . foldl' (\a b -> a `shiftL` 8 .|. fromIntegral b) 0 $
               take 8 $ (++ repeat 0) $
               B.zipWith (\b q -> shiftL (b .&. 0xE) 4 .|. (min 31 $ max 33 q - 33)) sq qs

fromTags :: String -> String -> BamRaw -> Index
fromTags itag qtag br = fromSQ sq  (if B.null qs then "@@@@@@@@" else qs)
  where
    sq = br_extAsString itag br
    qs = br_extAsString qtag br

gather :: MonadIO m => Refs -> Iteratee [BamRaw] m (U.Vector (Index, Index))
gather refs = filterStream (\b -> not (br_isPaired b) || br_isFirstMate b) =$
              progress "reading " (hPutStr stderr) refs =$
              mapStream (fromTags "XI" "YI" &&& fromTags "XJ" "YJ") =$
              subsam2vector 100000

subsam2vector :: (MonadIO m, ListLike s a, Nullable s, VG.Vector v a) => Int -> Iteratee s m (v a)
subsam2vector sz = liftIO (VGM.new sz) >>= go 0
  where
    go !i !mv = tryHead >>= \x -> case x of
                  Nothing -> liftIO $ if i < sz then VG.unsafeFreeze $ VGM.take i mv
                                                else VG.unsafeFreeze mv
                  Just  a -> do -- when (i `rem` 0x10000 == 0) $ liftIO performGC
                                liftIO $ if i < sz
                                    then VGM.write mv i a
                                    else do p <- randomRIO (0,i)
                                            when (p < sz) $ VGM.write mv p a
                                go (i+1) mv

data Both = Both { p7is :: U.Vector Index
                 , p7ns :: V.Vector T.Text
                 , p5is :: U.Vector Index
                 , p5ns :: V.Vector T.Text }

instance FromJSON Both where
    parseJSON = withObject "toplevel object expected" $ \v ->
                          both <$> ((v .: "p7index") >>= parse_assocs)
                               <*> ((v .: "p5index") >>= parse_assocs)
      where
        parse_assocs = withObject "association list expected" $ \o ->
                            sequence [ (,) k <$> withText "sequence expected" (return . T.encodeUtf8) v | (k,v) <- HM.toList o ]

        both as7 as5 = Both (U.fromList is7') (V.fromList ns7) (U.fromList is5') (V.fromList ns5)
          where
            (ns7,is7) = unzip as7
            (ns5,is5) = unzip as5
            is7' = map fromS is7
            is5' = map fromS is5


-- | Compute mismatch score: sum of the qualities in 'a' at positions
-- where the bases don't match.  Works by comparing through an xor,
-- building a mask from it, then adding quality scores sideways.
--
-- Since we keep quality scores in the lower 5 bits of each byte, adding
-- all eight is guaranteed to fit into the highest 8 bits.
match :: Index -> Index -> Word64
match (Index a) (Index b) = score
  where x = a `xor` b
        y = (shiftR x 5 .|. shiftR x 6 .|. shiftR x 7) .&. 0x0101010101010101
        mask = (0x2020202020202020 - y) .&. 0x1F1F1F1F1F1F1F1F
        score = shiftR ((a .&. mask) * 0x0101010101010101) 56

-- | A mixture description is one probability for each combination of p7
-- and p5 index.  They should sum to one.
type Mix = VS.Vector Double
type MMix = VSM.IOVector Double
padding = 31

-- | Computing the naively assumed mix when nothing is known:  uniform
-- distribution.
unknownMix :: (Int,Int) -> Mix
unknownMix (n7,n5) = VS.replicate total (recip $ fromIntegral total)
  where
    !total = n7*((n5+padding) .&. complement padding)

-- | Computing the naively assumed mix:  we assume 99% is composed of
-- an even mix of the provided list, the remaining 1% gets distributed
-- among everything else.
naiveMix :: (Int,Int) -> [(Int,Int)] -> Mix
naiveMix (n7,n5') ps = VS.replicate total (0.01 / fromIntegral (total-known))
                 VS.// [ (n5*i7+i5, 0.99 / fromIntegral known) | (i7,i5) <- ps ]
  where
    !n5 = (n5' + padding) .&. complement padding
    !total = n7*n5
    !known = length ps

-- Matches an index against both p7 and p5 lists, computes posterior
-- likelihoods from the provided prior.
unmix1 :: U.Vector Index -> U.Vector Index -> Mix -> MMix -> (Index, Index) -> IO ()
unmix1 p7 p5 prior acc (x,y) = do
    let !m7 = VS.fromListN (U.length p7) . map (phredPow . match x) $ U.toList p7
        !l5 = (U.length p5 + padding) .&. complement padding
        !m5 = VS.fromListN l5 $ map (phredPow . match y) (U.toList p5) ++ repeat 0

    -- *sigh*, concatMap doesn't fuse.  Gotta do it manually.  :(
    -- v <- VSM.new (VS.length m7 * l5)
    {- let loop !i !j !k !acc
            | i == VS.length m7 = loop2 (recip acc) 0
            | j == VS.length m5 = loop (i+1) 0 k acc
            | otherwise = let !p = prior VS.! k * m7 VS.! i * m5 VS.! j
                          in VSM.write v k p >> loop i (j+1) (k+1) (acc + p)

        loop2 !tot !k
            | k == VS.length m7 * VS.length m5 = return v
            | otherwise = VSM.read v k >>= VSM.write v k . (*) tot >> loop2 tot (k+1)

    loop 0 0 0 0 -}
    VSM.unsafeWith acc                  $ \pw ->
        VS.unsafeWith prior             $ \pv ->
            VS.unsafeWith m7            $ \p7 ->
                VS.unsafeWith m5        $ \p5 ->
                    c_loop pw pv p7 (VS.length m7) p5 (l5 `div` (1+padding))
    -- VS.unsafeFreeze v

foreign import ccall unsafe "c_loop"
    c_loop :: Ptr Double -> Ptr Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

phredPow :: Word64 -> Double
phredPow x = exp $ -0.1 * log 10 * fromIntegral x

-- | One iteration of the EM algorithm.  Input is a vector of pairs of
-- indices, the p7 and p5 index collections, and a prior mixture; output
-- is the posterior mixture.
iterEM :: U.Vector (Index, Index) -> U.Vector Index -> U.Vector Index -> Mix -> IO Mix
iterEM pairs p7 p5 prior = do
    acc <- VSM.replicate (VS.length prior) 0
    U.mapM_ (unmix1 p7 p5 prior acc) pairs
    VS.unsafeFreeze acc

main :: IO ()
main = do
    Just Both{..} <- fmap decodeStrict' $ B.readFile . (++ "/usr/share/lims/global_index_list.json") =<< getHomeDirectory
    print $ U.length p7is
    print $ U.length p5is

    ixvec <- concatDefaultInputs >=> run $ gather . meta_refs
    print $ U.length ixvec

    let loop v = do v' <- iterEM ixvec p7is p5is v
                    inspect p7ns p5ns v'
                    loop v'
    loop $ unknownMix (U.length p7is, U.length p5is)

inspect :: V.Vector T.Text -> V.Vector T.Text -> Mix -> IO ()
inspect n7 n5 mix = do
    getClockTime >>= print
    v <- U.unsafeThaw $ U.fromListN (VS.length mix) $ zip [0..] $ VS.toList mix
    V.partialSortBy (\(_,a) (_,b) -> compare b a) v 20
    v' <- U.unsafeFreeze v
    U.forM_ (U.take 20 v') $ \(i,n) -> do
       let (i7, i5) = i `quotRem` ((V.length n5 + padding) .&. complement padding)
       T.print "{}, {}: {}\n"
            ( T.left 7 ' ' $ n7 V.! i7
            , T.left 7 ' ' $ n5 V.! i5
            , T.left 8 ' ' $ T.fixed 2 n )
    putStrLn []

