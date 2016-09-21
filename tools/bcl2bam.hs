{-# LANGUAGE RecordWildCards, PatternGuards #-}
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


-- XXX  Single read, for now.  It needs the read length definitions to
-- be correct.
tileToBam :: String -> Tile -> [BamRec]
tileToBam expnm Tile{ tile_locs = Locs vlocs, ..} =
    V.ifoldr one_cluster [] vlocs

  where
    one_cluster i (px,py) = (:) nullBamRec {
        b_qname = fromString (printf "%s:%d:%d:%d" expnm tile_nbr px py),
        b_flag  = flagUnmapped,
        b_seq   = V.convert $ V.map (get_nucs i) tile_bcls,
        b_qual  = V.convert $ V.map (get_qual i) tile_bcls,
        b_exts  = [] }

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
          else let cycles_read_one = (minimum cycles, maximum cycles)
               in return $ Just LaneDef{..}
      else return Nothing

  where
    -- we deal with every tile we can find a locs file for
    get_tile fn | takeExtension fn == ".gz"                  = get_tile (dropExtension fn)
                | takeExtension fn `elem` [".clocs",".locs"] = Just $ dropExtension fn
                | otherwise                                  = Nothing

    get_cycle ('C':fn) | (l,".1") <- break (== '.') fn, [(c,[])] <- reads l = Just c
    get_cycle        _                                                      = Nothing



-- For each tile, read locs and all the bcls.  Run tileToBam and emit.
bamFromBcl :: LaneDef -> Enumerator [BamRec] IO b
bamFromBcl LaneDef{..} iter0 =
    foldM (\iter fn -> do tile <- one_tile fn
                          enumPure1Chunk (tileToBam experiment tile) iter)
          iter0 tiles

  where
    (c1,ce) = cycles_read_one

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

        get_bcls = fmap V.fromList . forM [c1..ce] $ \ncycle -> do
            let fn_bcl = path_bcl ++ "/C" ++ show ncycle ++ ".1/" ++ fn ++ ".bcl"
            try_read_or (fn_bcl <.> "gz") readBCL $
                try_read_or fn_bcl readBCL $
                    return (BCL V.empty)

        try_read_or f r k = do e <- doesFileExist f
                               if e then r f else k


main = do
    lanedefs <- getArgs >>= mapM lanesFromRun >>= return . concat
    foldr ((>=>) . bamFromBcl) run lanedefs $ pipeBamOutput mempty
