{-# LANGUAGE RecordWildCards, PatternGuards, Arrows #-}
-- BCL to BAM command line driver
-- This should become part of Jivebunny; we'll see if the standalone
-- version will retain any value.
--
-- Basic idea:  Input is a directory of BCL files, a directory of
-- pos/locs/clocs files, read length definitions.  Output is BAM.  (An
-- argument could be made that we also need tiles, swathes and what
-- not.)
--
-- If instead we start from a run folder, we can get the read length
-- definitions from the RunInfo.xml file, and the directories are
-- implied.  We need to specify a lane in this case.
--
-- In practice, we might want to start from the run folder and allow
-- everything to be overridden, or override everything from a clean
-- slate.


import Bio.Bam
import Bio.Prelude
import System.Directory
import System.FilePath
import Text.XML.HXT.Core

import qualified Data.ByteString as B
import qualified Data.Set as S
import qualified Data.Vector as VV
import qualified Data.Vector.Generic as V

import Bio.Illumina.BCL
import Bio.Illumina.Locs

-- conversion of BCL/LOCS to BAM.  We go tile by tile, and for each tile
-- we need a bunch of BCLs (one per cycle) and a LOCS.  (Note on
-- parallelization:  we can read the many files in parallel, and we
-- probably should.  We can also read the files for one tile while
-- outputting another.)

data Tile = Tile
    { tile_nbr :: !Int
    , tile_locs :: !Locs
    , tile_bcls :: !(VV.Vector BCL) }


tileToBam :: LaneDef -> Tile -> [[BamRec]]
tileToBam LaneDef{..} Tile{ tile_locs = Locs vlocs, ..}
    = zipWith one_cluster [0..] (V.toList vlocs)
  where
    -- XXX  Need to generate up to two index reads and up to two BAM records.
    one_cluster i (px,py) =
        nullBamRec { b_qname = qname,
                     b_flag  = maybe flagsSingle (const flagsReadOne) cycles_read_two,
                     b_seq   = V.convert $ get_seq cycles_read_one,
                     b_qual  = V.convert $ get_quals cycles_read_one,
                     b_exts  = common_exts } :
        case cycles_read_two of
            Nothing -> []
            Just r2 -> nullBamRec { b_qname = qname,
                                    b_flag  = flagsReadTwo,
                                    b_seq   = V.convert $ get_seq r2,
                                    b_qual  = V.convert $ get_quals r2,
                                    b_exts  = common_exts } : []
      where
        qname       = fromString (printf "%s:%d:%d:%d" experiment tile_nbr px py)
        common_exts = maybe id (indexRead "XI" "YI") cycles_index_one $
                      maybe id (indexRead "XJ" "YJ") cycles_index_two $ []

        get_seq   (ra,re) = V.map (get_nucs i) $ V.slice (ra-1) (re-ra+1) tile_bcls
        get_quals (ra,re) = V.map (get_qual i) $ V.slice (ra-1) (re-ra+1) tile_bcls

        indexRead k1 k2 rng
            = insertE k1 (Text . fromString . V.toList . V.map showNucleotides          $ get_seq   rng)
            . insertE k2 (Text . B.pack . V.toList . V.map ((+33) . fromIntegral . unQ) $ get_quals rng)

    flagsSingle = flagUnmapped
    flagsReadOne = flagUnmapped .|. flagMateUnmapped .|. flagPaired .|. flagFirstMate
    flagsReadTwo = flagUnmapped .|. flagMateUnmapped .|. flagPaired .|. flagSecondMate

    get_qual i (BCL v) = Q . maybe 0 (`shiftR` 2) $ v V.!? i
    get_nucs i (BCL v) = maybe nucsN code_to_nucs $ v V.!? i
    code_to_nucs x |    x    == 0 = nucsN
                   | x .&. 3 == 0 = nucsA
                   | x .&. 3 == 1 = nucsC
                   | x .&. 3 == 2 = nucsG
                   | otherwise    = nucsT


-- | Definition of a lane to be processed.  Includes paths, read
-- definitions.
data LaneDef = LaneDef {
    experiment :: String,
    lane_number :: Int,

    -- | Root of BCL hierarchy, contains BCLs in subdirectories, filter
    -- and control files.
    path_bcl :: FilePath,

    -- | Path to location files.  Contains clocs, locs, or pos_txt.
    path_locs :: FilePath,

    -- | Cycles in the first business read.
    cycles_read_one :: (Int,Int),

    -- | Cycles in the second business read, if present.
    cycles_read_two :: Maybe (Int,Int),

    -- | Cycles in the first index read, if present.
    cycles_index_one :: Maybe (Int,Int),

    -- | Cycles in the second index read, if present.
    cycles_index_two :: Maybe (Int,Int),

    -- | List of basenames, one for each tile.
    tiles :: [String] }
  deriving Show

-- | Takes a run folder and derives the lane definitions.
lanesFromRun :: FilePath -> IO [LaneDef]
lanesFromRun rundir = fmap catMaybes . forM [1..8] $ \lane_number -> do
    let experiment = case break (== '_') (takeBaseName rundir) of       -- try and drop the date
                        (l,r) | all isDigit l -> drop 1 r
                        _                     -> takeBaseName rundir

    let path_bcl = rundir </> "Data/Intensities/BaseCalls/L00" ++ [intToDigit lane_number]
        path_locs = rundir </> "Data/Intensities/L00" ++ [intToDigit lane_number]

    has_both <- (&&) <$> doesDirectoryExist path_bcl
                     <*> doesDirectoryExist path_locs

    if has_both then do
        tiles <- S.toList . S.fromList . mapMaybe get_tile <$> getDirectoryContents path_locs
        cycles <- mapMaybe get_cycle <$> getDirectoryContents path_bcl
        if null cycles || null tiles
          then return Nothing
          else do ri <- expand_ri (maximum cycles) 1 <$> getRunInfo rundir
                  let (cycles_read_one, cycles_read_two)
                        = case map fst $ filter ((Just True /=) . snd) ri of
                            r1:r2:_ -> (r1, Just r2)
                            r1:_    -> (r1, Nothing)
                            _       -> error "shouldn't happen"

                  let (cycles_index_one, cycles_index_two)
                        = case map fst $ filter ((Just False /=) . snd) ri of
                            r1:r2:_ -> (Just r1, Just r2)
                            r1:_    -> (Just r1, Nothing)
                            _       -> (Nothing, Nothing)

                  return $ Just LaneDef{..}

      else return Nothing

  where
    -- we deal with every tile we can find a locs file for
    get_tile fn | takeExtension fn == ".gz"                  = get_tile (dropExtension fn)
                | takeExtension fn `elem` [".clocs",".locs"] = Just $ dropExtension fn
                | otherwise                                  = Nothing

    get_cycle ('C':fn) | (l,".1") <- break (== '.') fn, [(c,[])] <- reads l = Just c
    get_cycle        _                                                      = Nothing

    expand_ri total count [           ] = if count <= total then [((count,total),Nothing)] else []
    expand_ri total count ((l,isix):rs) = ((count,count+l-1),Just isix) : expand_ri total (count+l) rs


-- For each tile, read locs and all the bcls.  Run tileToBam and emit.
bamFromBcl :: LaneDef -> Enumerator [BamRec] IO b
bamFromBcl ld@LaneDef{..} iter0 =
    foldM (\iter fn -> do tile <- one_tile fn
                          enumList (tileToBam ld tile) iter)
          iter0 tiles
  where
    !ce = maybe id (max . snd) cycles_index_two $
          maybe id (max . snd) cycles_index_one $
          maybe id (max . snd) cycles_read_two  $
          snd cycles_read_one

    one_tile :: FilePath -> IO Tile
    one_tile fn = Tile nbr <$> get_locs <*> get_bcls
      where
        nbr = case reads . reverse . takeWhile (/= '_') . reverse $ fn of
                    [(n,"")] -> n ; _ -> 0

        get_locs = do
            let fn_locs = path_locs </> fn
            try_read_or (fn_locs <.> ".clocs.gz") readClocs $
                try_read_or (fn_locs <.> ".clocs") readClocs $
                    try_read_or (fn_locs <.> ".locs.gz") readLocs $
                        try_read_or (fn_locs <.> ".locs") readLocs $
                            return (Locs V.empty)

        get_bcls = fmap V.fromList . forM [1..ce] $ \ncycle -> do
            let fn_bcl = path_bcl ++ "/C" ++ show ncycle ++ ".1/" ++ fn ++ ".bcl"
            try_read_or (fn_bcl <.> "gz") readBCL $
                try_read_or fn_bcl readBCL $
                    return (BCL V.empty)

        try_read_or f r k = do e <- doesFileExist f
                               if e then r f else k


main  = do
    rs <- getArgs

    lanedefs <- mapM lanesFromRun rs >>= return . concat
    hPutStrLn stderr $ show lanedefs
    foldr ((>=>) . bamFromBcl) run lanedefs $ pipeBamOutput mempty


-- Look for a useable XML file, either RunInfo.xml, or
-- RunParameters.xml.  We'll match case insensitively, because
-- sometimes case gets mangled during the network copy.
getRunInfo :: FilePath -> IO [(Int,Bool)]
getRunInfo dir = do
    xmls <- filter (\f -> map toLower f == "runinfo.xml" || map toLower f == "runparameters.xml")
            <$> getDirectoryContents dir
    case xmls of
        fp:_ -> fmap (map snd . sort) .
                runX $ readDocument [ withValidate no ] (dir </> fp)
                       >>> deep (isElem >>> hasName "Reads")
                       >>> getChildren >>> hasName "Read"
                       >>> toReadDef ()
        [  ] -> return []

toReadDef :: ArrowXml a => t -> a XmlTree (Int, (Int, Bool))
toReadDef _ = ( getAttrValue "Number" >>> readA ) &&&
              ( getAttrValue "NumCycles" >>> readA ) &&&
              ( getAttrValue "IsIndexedRead" >>> arr boolF )
  where
    readA = arrL $ \s -> case reads (dropWhile isSpace s) of
                            [(a,b)] | all isSpace b -> [a] ; _ -> []
    boolF s = if s == "Y" || s == "y" then True else False
