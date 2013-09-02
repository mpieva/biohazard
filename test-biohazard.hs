import Test.HUnit
import Bio.File.Bgzf

import Bio.Iteratee
import System.Posix.Types ( COff )
import System.Exit
import qualified Data.ByteString as S

main :: IO ()
main = do cnts <- runTestTT $ TestList [ TestLabel "BGZF" test_bgzf ]
          if failures cnts == 0 then exitSuccess else exitFailure

test_bgzf :: Test
test_bgzf = TestList [ TestCase check_bam_block_list, TestCase check_bam_block_list_seek ]

check_bam_block_list :: Assertion
check_bam_block_list =
    -- fileDriver (joinI $ decompress' (virtualSeek 0 >> print_block))
    fileDriver (joinI $ decompress' $ getBlockList ) some_bgzf_file
    >>= assertEqual "unexpected block list" known_blocks

check_bam_block_list_seek :: Assertion
check_bam_block_list_seek =
    fileDriverRandom (joinI $ decompress' $ virtualSeek 0 >> getBlockList ) some_bgzf_file
    >>= assertEqual "unexpected block list" known_blocks


some_bgzf_file :: FilePath
some_bgzf_file = "/mnt/ngs_data/101203_SOLEXA-GA04_00007_PEDi_MM_QF_SR/Ibis/BWA/s_5_L3280_sequence_mq_hg19_nohap.bam"

known_blocks :: [(COff, Int)]
known_blocks = [(0,65536), (1627521024,65536), (3209297920,65536), (4116643840,65536), (5036769280,65536), (6071189504,65536),
    (7127957504,65536), (8169455616,65536), (9189064704,65536), (10195501056,65536), (11269373952,65536), (12313034752,65536),
    (13372882944,65536), (14407958528,65536), (15697903616,65536), (17157521408,65536), (18614190080,65536), (20095696896,65536),
    (21540044800,65536), (23005495296,65536), (24441257984,65536), (25919488000,65536), (27378712576,65536), (28806414336,65536),
    (30298144768,65536), (31775916032,65536), (33248444416,65536), (34703212544,65536), (36106272768,65536), (37567791104,65536),
    (39027146752,65536), (40473722880,65536), (41894150144,65536), (43395317760,65536), (44840124416,65536), (46291025920,65536),
    (47691923456,65536), (49136533504,65536), (50582781952,65536), (52029947904,65536), (53468528640,65536), (54865559552,65536),
    (56325046272,65536), (57735839744,65536), (59189428224,65536), (60636659712,65536), (62072291328,65536), (63498092544,65536),
    (64939360256,65536), (66396291072,65536), (67904864256,65536), (69389516800,65536), (70850248704,65536), (72320876544,65536),
    (73770729472,65536), (75232051200,65536), (76726927360,65536), (78167736320,65536), (79583838208,65536), (80989716480,65536),
    (82401034240,39996), (83302875136,0)]

getBlockList :: Iteratee Block IO [(COff, Int)]
getBlockList = liftI $ step []
  where
    step acc (Chunk (Block p s)) = let !l = S.length s in liftI $ step ((p,l):acc)
    step acc e@(EOF _) = idone (reverse acc) e


