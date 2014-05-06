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


{- moved from Codec.BGZF and Bio.Bam.Rec

import Bio.Bam
import Codec.Bgzf
import System.IO

import Data.ByteString ( hPut )
import qualified Data.ByteString.Char8 as S

some_file :: FilePath
some_file = "/mnt/ngs_data/101203_SOLEXA-GA04_00007_PEDi_MM_QF_SR/Ibis/BWA/s_5_L3280_sequence_mq_hg19_nohap.bam"

bgzf_test :: FilePath -> IO ()
bgzf_test = fileDriver $
            joinI $ decompressBgzfBlocks $
            mapChunksM_ (print . block_offset)

bam_test' :: FilePath -> IO ()
bam_test' = fileDriver $
            joinI $ decompressBgzfBlocks $
            joinI $ decodeBam      $
            dump_bam

bam_test :: FilePath -> IO ()
bam_test = fileDriverRandom $
           joinI $ decompressBgzfBlocks $
           joinI $ do -- decodeBam $ \h -> joinI $ takeStream 100 (dump_bam h)
                      seek 0
                      decodeBam dump_bam

dump_bam :: BamMeta -> Iteratee [BamRaw] IO ()
dump_bam meta = lift (print meta) >> print_names

seek_test :: [Char] -> Word32 -> IO ()
seek_test fp i = do
    idx <- readBamIndex $ fp ++ ".bai"
    flip fileDriverRandom fp $
           joinI $ decompressBgzfBlocks $
           joinI $ decodeBamSequence idx (Refseq i) print_names_and_refs

sam_test :: IO ()
sam_test = fileDriver (joinI $ decodeSam (const print_names')) "foo.sam"

print_names :: Iteratee [BamRaw] IO ()
print_names = mapStreamM_ $ S.putStrLn . br_qname

print_names_and_refs :: Iteratee [BamRaw] IO ()
print_names_and_refs = mapStreamM_ pr
  where pr b = putStrLn $ shows (br_qname b) " " ++ show (br_rname b)

print_names' :: Iteratee [BamRec] IO ()
print_names' = mapStreamM_ $ S.putStrLn . b_qname


bam2bam_test :: IO ()
bam2bam_test = withFile "foo.bam" WriteMode $       \hdl ->
               flip fileDriver some_file    $
               joinI $ decompressBgzfBlocks $
               joinI $ decodeBam            $       \meta ->
               joinI $ encodeBam meta       $
               mapChunksM_ (hPut hdl)

sam2bam_test :: IO ()
sam2bam_test = withFile "bar.bam" WriteMode       $             \hdl ->
               flip fileDriver "foo.sam"          $
               joinI $ decodeSam                  $             \meta ->
               joinI $ mapStream encodeBamEntry   $
               lift (print meta)                >>=             \_ ->
               joinI $ encodeBam meta             $
               mapChunksM_ (hPut hdl)
-}


