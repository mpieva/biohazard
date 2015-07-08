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
--  - Input layer to gather index sequences.  (Done.)
--  - Input layer to gather read group defs.  (Done.)
--  - First pass to gather data.  Any index read shall be represented
--    in a single Word64.  (Done.  Reading BAM is slow.)
--  - Multiple passes of the EM algorithm.  (Done.)
--  - Start with a naive mix, to avoid arguments.  (Done.)
--  - Final calling pass from BAM to BAM.  (Done.)
--  - Auxillary statistics:  composition of the mix, false assignment
--    matrix or rates per read group, maximum achievable false
--    assignment rates.

import Bio.Bam
import Bio.Util ( showNum )
import Control.Arrow ( (&&&) )
import Control.Monad ( when, unless, forM_ )
import Data.Aeson
import Data.Bits
import Data.Char ( chr )
import Data.List ( foldl' )
import Data.Monoid
import Data.Vector.Unboxed.Deriving
import Data.Version ( showVersion )
import Data.Word ( Word64 )
import Foreign.C.Types
import Foreign.Marshal.Alloc
import Foreign.Ptr
import Foreign.Storable
import Paths_biohazard ( version )
import System.Console.GetOpt
import System.Directory ( getHomeDirectory )
import System.Environment ( getProgName, getArgs )
import System.Exit
import System.IO
import System.IO.Unsafe ( unsafePerformIO )
import System.Random ( randomRIO )
import System.Time ( getClockTime )

import qualified Data.ByteString as B
import qualified Data.HashMap.Strict as HM
import qualified Data.Map as M
import qualified Data.Text as T
import qualified Data.Text.Encoding as T
import qualified Data.Text.Format as T
import qualified Data.Text.IO as T
import qualified Data.Text.Lazy.Builder as TB
import qualified Data.Text.Lazy.Builder.Int as TB
import qualified Data.Vector as V
import qualified Data.Vector.Algorithms.Intro as V
import qualified Data.Vector.Unboxed as U
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

fromS :: B.ByteString -> Index
fromS sq = fromSQ sq (B.replicate (B.length sq) 64)

fromSQ :: B.ByteString -> B.ByteString -> Index
fromSQ sq qs = Index . foldl' (\a b -> a `shiftL` 8 .|. fromIntegral b) 0 $
               take 8 $ (++ repeat 0) $
               B.zipWith (\b q -> shiftL (b .&. 0xE) 4 .|. (min 31 $ max 33 q - 33)) sq qs

fromTags :: String -> String -> BamRaw -> Index
fromTags itag qtag br = fromSQ sq  (if B.null qs then "@@@@@@@@" else qs)
  where
    sq = br_extAsString itag br
    qs = br_extAsString qtag br

gather :: MonadIO m => Int -> (String -> IO ()) -> (String -> IO ()) -> BamMeta -> Iteratee [BamRaw] m (U.Vector (Index, Index))
gather num say mumble hdr = case hdr_sorting $ meta_hdr hdr of
    Unsorted    -> greedy
    Grouped     -> greedy
    Queryname   -> greedy
    Unknown     -> safe
    Coordinate  -> fair
    GroupSorted -> fair
  where
    greedy = do liftIO . say $ "File is unsorted, sampling up to "
                            ++ showNum num ++ " records from the beginning.\n"
                go stream2vectorN

    fair   = do liftIO . say $ "File is sorted, need to sample up to "
                            ++ showNum num ++ " from whole file.\n"
                go subsam2vector

    safe   = do liftIO . say $ "File might be sorted, need to sample up to "
                            ++ showNum num ++ " from whole file.\n"
                go subsam2vector

    go k = filterStream (\b -> not (br_isPaired b) || br_isFirstMate b) =$
           progress "reading " mumble (meta_refs hdr) =$
           mapStream (fromTags "XI" "YI" &&& fromTags "XJ" "YJ") =$ k num


subsam2vector :: (MonadIO m, ListLike s a, Nullable s, VG.Vector v a) => Int -> Iteratee s m (v a)
subsam2vector sz = liftIO (VGM.new sz) >>= go 0
  where
    go !i !mv = tryHead >>= \x -> case x of
                  Nothing -> liftIO $ if i < sz then VG.unsafeFreeze $ VGM.take i mv
                                                else VG.unsafeFreeze mv
                  Just  a -> do liftIO $ if i < sz
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

data RG = RG { rgid :: B.ByteString
             , rgi7 :: Int
             , rgi5 :: Int
             , tags :: BamOtherShit }

-- | Parses read group defintions from a file.  The file can have
-- optional header lines, the remainder must be a tab-separated table,
-- first column is the read group name, second is the P7 index name,
-- third is the P5 index name (*must* be present), all others are tagged
-- fields just like BAM expects them in the header.
--
-- For integration with a LIMS, something structured like JSON would
-- probably work better, however, absent such a LIMS, tables are easier
-- to come by.

readRGdefns :: V.Vector T.Text -> V.Vector T.Text -> T.Text -> [ RG ]
readRGdefns p7ns p5ns = map repack . filter (not . null) . map (T.split (=='\t'))
                      . dropWhile ("#" `T.isPrefixOf`) . T.lines
  where
    repack (rg:p7:p5:tags) = case V.elemIndex p7 p7ns of
        Nothing -> error $ "unknown P7 index " ++ show p7
        Just i7 -> case V.elemIndex p5 p5ns of
            Nothing -> error $ "unknown P5 index " ++ show p5
            Just i5 -> RG (T.encodeUtf8 rg) i7 i5 (map repack1 tags)
    repack ws = error $ "short RG line " ++ show (T.intercalate "\t" ws)
    repack1 w | T.length w > 3 && T.index w 2 == ':' = (T.index w 1, T.index w 2, T.encodeUtf8 $ T.drop 3 w)
              | otherwise = error $ "illegal tag " ++ show w

default_rgs :: T.Text
default_rgs = "PhiXA\tPhiA\tPhiA\nPhiXC\tPhiC\tPhiC\nPhiXG\tPhiG\tPhiG\nPhiXT\tPhiT\tPhiT\nPhiX\tPhiX\tis4\n"

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

padding :: Int
padding = 31

-- | Computing the naively assumed mix when nothing is known:  uniform
-- distribution.
naiveMix :: (Int,Int) -> Int -> Mix
naiveMix (n7,n5) total = VS.replicate vecsize (fromIntegral total / fromIntegral bins)
  where
    !vecsize = n7 * ((n5+padding) .&. complement padding)
    !bins    = n7 * n5

-- | Matches an index against both p7 and p5 lists, computes posterior
-- likelihoods from the provided prior and accumulates them onto the
-- posterior.
unmix1 :: U.Vector Index -> U.Vector Index -> Mix -> MMix -> (Index, Index) -> IO ()
unmix1 p7 p5 prior acc (x,y) = do
    let !m7 = VS.fromListN (U.length p7) . map (phredPow . match x) $ U.toList p7
        !l5 = (U.length p5 + padding) .&. complement padding
        !m5 = VS.fromListN l5 $ map (phredPow . match y) (U.toList p5) ++ repeat 0

    -- *sigh*, Vector doesn't fuse well.  Gotta hand it over to gcc.  :-(
    VSM.unsafeWith acc                  $ \pw ->
        VS.unsafeWith prior             $ \pv ->
            VS.unsafeWith m7            $ \q7 ->
                VS.unsafeWith m5        $ \q5 ->
                    c_unmix_acc1 pw pv q7 (fromIntegral $ VS.length m7)
                                       q5 (fromIntegral $ l5 `div` succ padding)

foreign import ccall unsafe "c_unmix_acc1"
    c_unmix_acc1 :: Ptr Double ->               -- accumulator (mutable)
                    Ptr Double ->               -- prior mix
                    Ptr Double -> CInt ->       -- P7 scores, length
                    Ptr Double -> CInt -> IO () -- P5 scores, length

-- | Matches an index against both p7 and p5 lists, computes MAP
-- assignment and quality score.
class1 :: U.Vector Index -> U.Vector Index -> Mix -> (Index, Index) -> (Double,Int,Int)
class1 p7 p5 prior (x,y) =
    let !m7 = VS.fromListN (U.length p7) . map (phredPow . match x) $ U.toList p7
        !l5 = (U.length p5 + padding) .&. complement padding
        !m5 = VS.fromListN l5 $ map (phredPow . match y) (U.toList p5) ++ repeat 0

    -- *sigh*, Vector doesn't fuse well.  Gotta hand it over to gcc.  :-(
    in unsafePerformIO              $
       alloca                       $ \pi7 ->
       alloca                       $ \pi5 ->
       VS.unsafeWith prior          $ \pv ->
       VS.unsafeWith m7             $ \q7 ->
       VS.unsafeWith m5             $ \q5 ->
       (,,) <$> c_class1 pv q7 (fromIntegral $ VS.length m7)
                            q5 (fromIntegral $ l5 `div` succ padding) pi7 pi5
            <*> (fromIntegral <$> peek pi7) <*> (fromIntegral <$> peek pi5)

foreign import ccall unsafe "c_class1"
    c_class1 :: Ptr Double ->               -- prior mix
                Ptr Double -> CInt ->       -- P7 scores, length
                Ptr Double -> CInt ->       -- P5 scores, length
                Ptr CInt -> Ptr CInt ->     -- OUT i7, i5
                IO Double                   -- posterior prob(?)

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

data Loudness = Quiet | Normal | Loud

data Conf = Conf {
        cf_index_list :: FilePath,
        cf_output     :: Maybe (BamMeta -> Iteratee [BamRaw] IO ()),
        cf_threshold  :: Double,
        cf_loudness   :: Loudness,
        cf_readgroups :: [FilePath] }

defaultConf :: IO Conf
defaultConf = do home <- getHomeDirectory
                 return $ Conf {
                        cf_index_list = home ++ "/usr/share/lims/global_index_list.json",
                        cf_output     = Nothing,
                        cf_threshold  = 0.000005,
                        cf_loudness   = Normal,
                        cf_readgroups = [] }

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o" ["output"]         (ReqArg set_output   "FILE") "Send output to FILE",
    Option "I" ["index-database"] (ReqArg set_index_db "FILE") "Read index database from FILE",
    Option "r" ["read-groups"]    (ReqArg set_rgs      "FILE") "Read read group definitions from FILE",
    Option [ ] ["threshold"]      (ReqArg set_thresh   "FRAC") "Iterate till uncertainty is below FRAC",
    Option "v" ["verbose"]        (NoArg             set_loud) "Print more diagnostic messages",
    Option "q" ["quiet"]          (NoArg            set_quiet) "Print fewer diagnostic messages",
    Option "h?" ["help", "usage"] (NoArg        (const usage)) "Print this message and exit",
    Option "V"  ["version"]       (NoArg         (const vrsn)) "Display version number and exit" ]
  where
    set_output  "-" c = return $ c { cf_output = Just $ pipeRawBamOutput }
    set_output   fp c = return $ c { cf_output = Just $ writeRawBamFile fp }
    set_index_db fp c = return $ c { cf_index_list = fp }
    set_rgs      fp c = return $ c { cf_readgroups = fp : cf_readgroups c }
    set_loud        c = return $ c { cf_loudness = Loud }
    set_quiet       c = return $ c { cf_loudness = Quiet }
    set_thresh    a c = readIO a >>= \x -> return $ c { cf_threshold = x }

    usage = do pn <- getProgName
               putStrLn $ usageInfo ("Usage: " ++ pn ++ " [options] [bam-files]\n" ++
                                     "Decomposes a mix of libraries and assigns read groups.") options
               exitSuccess

    vrsn = do pn <- getProgName
              hPutStrLn stderr $ pn ++ ", version " ++ showVersion version
              exitSuccess



main :: IO ()
main = do
    (opts, files, errs) <- getOpt Permute options <$> getArgs
    unless (null errs) $ mapM_ (hPutStrLn stderr) errs >> exitFailure
    when (null files) $ hPutStrLn stderr "no input files." >> exitFailure
    Conf{..} <- foldl (>>=) defaultConf opts
    add_pg <- addPG $ Just version

    let notice  = case cf_loudness of Quiet -> \_ -> return () ; _ -> hPutStr stderr
        info    = case cf_loudness of Loud  -> hPutStr stderr ;  _ -> \_ -> return ()

    Just Both{..} <- fmap decodeStrict' $ B.readFile cf_index_list
    rgdefs <- concatMap (readRGdefns p7ns p5ns) . (:) default_rgs <$> mapM T.readFile cf_readgroups
    notice $ "Got " ++ showNum (U.length p7is) ++ " P7 indices and "
                    ++ showNum (U.length p5is) ++ " P5 indices.\n"
    notice $ "Declared " ++ showNum (length rgdefs) ++ " read groups.\n"

    let !rgs     = HM.fromList [ ((i7,i5),rg) | RG !rg !i7 !i5 _ <- rgdefs ]
        inspect = inspect' rgs p7ns p5ns

    ixvec <- concatInputs files >=> run $ gather 50000 notice info
    notice $ "Got " ++ showNum (U.length ixvec) ++ " index pairs.\n"

    notice "decomposing mix "
    let loop !n v = do v' <- iterEM ixvec p7is p5is v
                       case cf_loudness of Loud   -> inspect v'
                                           Normal -> hPutStr stderr "."
                                           Quiet  -> return ()
                       let d = VS.foldl' (\a -> max a . abs) 0 $ VS.zipWith (-) v v'
                       if n > 0 && d > cf_threshold * fromIntegral (U.length ixvec)
                            then loop (n-1) v'
                            else do notice (if n == 0 then "\nmaximum number of iterations reached.\n"
                                                      else "\nmixture ratios converged.\n")
                                    return v'

    mix <- loop (50::Int) $ naiveMix (U.length p7is, U.length p5is) (U.length ixvec)

    case cf_loudness of
        Quiet  -> return ()
        _ -> do inspect mix
                hPutStrLn stderr "\nmaximum achievable scores:"
                forM_ rgdefs $ \RG{..} ->
                    let (p,_,_) = class1 p7is p5is mix (p7is U.! rgi7, p5is U.! rgi5)
                        q = negate . round $ 10 / log 10 * log p :: Int
                    in T.hprint stderr "{}: {}\n"
                            ( T.left 12 ' ' $ T.decodeUtf8 rgid
                            , T.left 3 ' ' $ TB.singleton 'Q' <> TB.decimal q )

    case cf_output of
        Nothing  -> return ()
        Just out -> concatInputs files >=> run $ \hdr ->
                        let hdr' = hdr { meta_other_shit =
                                          [ os | os@(x,y,_) <- meta_other_shit hdr, x /= 'R' || y /= 'G' ] ++
                                          HM.elems (HM.fromList [ (rgid, ('R','G', ('I','D',rgid):tags)) | RG{..} <- rgdefs ] ) }
                        in mapStream (\br ->
                                let x = fromTags "XI" "YI" br
                                    y = fromTags "XJ" "YJ" br
                                    (p,i7,i5) = class1 p7is p5is mix (x,y)
                                    q = negate . round $ 10 / log 10 * log p
                                    b = decodeBamEntry br
                                    b' = b { b_exts = M.delete "Z0" . M.delete "Z2" . M.insert "Z1" (Int q)
                                                    $ case HM.lookup (i7,i5) rgs of
                                                        Nothing -> M.delete "RG" $ b_exts b
                                                        Just rg -> M.insert "RG" (Text rg) $ b_exts b }
                                in encodeBamEntry b') =$
                           progress "writing " info (meta_refs hdr) =$
                           out (add_pg hdr')


inspect' :: HM.HashMap (Int,Int) B.ByteString -> V.Vector T.Text -> V.Vector T.Text -> Mix -> IO ()
inspect' rgs n7 n5 mix = do
    hPutStrLn stderr []
    getClockTime >>= hPutStrLn stderr . show
    v <- U.unsafeThaw $ U.fromListN (VS.length mix) $ zip [0..] $ VS.toList mix
    V.partialSortBy (\(_,a) (_,b) -> compare b a) v 20
    v' <- U.unsafeFreeze v
    U.forM_ (U.take 20 v') $ \(i,n) -> do
       let (i7, i5) = i `quotRem` ((V.length n5 + padding) .&. complement padding)
       T.hprint stderr "{}, {}: {} {} \n"
            ( T.left 7 ' ' $ n7 V.! i7
            , T.left 7 ' ' $ n5 V.! i5
            , T.left 8 ' ' $ T.fixed 2 n
            , case HM.lookup (i7,i5) rgs of
                Nothing -> ""
                Just rg -> T.concat [ "(", T.decodeUtf8 rg, ")" ] )

