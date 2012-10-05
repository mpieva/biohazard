import Control.Monad                        ( unless, foldM, (>=>) )
import Data.Iteratee.Base
import Data.Iteratee.IO
import Data.Iteratee.Iteratee
import Data.Iteratee.ListLike               ( mapStream )
import Data.Word                            ( Word8 )
import Bio.File.Bam
import Bio.File.Bam.Trim
import Bio.File.Bgzf
import Bio.Iteratee
import qualified Data.ByteString      as S  ( hPut )
import System.Console.GetOpt
import System.Environment                   ( getArgs )
import System.Exit                          ( exitFailure, exitSuccess )
import System.IO                            ( stdin, stdout, stderr, hPutStrLn )

data Conf = Conf { c_trim_pred :: [Nucleotide] -> [Word8] -> Bool
                 , c_pass_pred :: BamRec -> Bool }

options :: [OptDescr (Conf -> IO Conf)]
options = [ Option "q" ["minq"]   (ReqArg set_minq "Q") "Trim where quality is below Q"
          , Option "m" ["mapped"] (NoArg set_monly)     "Trim only mapped sequences"
          , Option "h?" ["help"]  (NoArg usage)         "Display this text" ]

set_minq :: String -> Conf -> IO Conf
set_minq s c = readIO s >>= \q -> return $ c { c_trim_pred = trim_low_quality q }

set_monly :: Conf -> IO Conf
set_monly c = return $ c { c_pass_pred = \r -> isMerged r || isUnmapped r }

usage :: Conf -> IO Conf
usage _ = do hPutStrLn stderr $ usageInfo info options ; exitSuccess
  where info = "Simple trimming of sequences in Bam files.  Reads a Bam file from stdin,\n\
               \trims sequences of low quality, writes Bam to stdout.  Does not trim\n\
               \merged reads."


main :: IO ()
main = do
    (opts, files, errors) <- getOpt Permute options `fmap` getArgs

    unless (null errors) $ mapM_ (hPutStrLn stderr) errors
    c <- foldM (flip id) (Conf (trim_low_quality 20) isMerged) opts
    unless (null errors && null files) exitFailure

    let do_trim r | c_pass_pred c r' = raw_data r
                  | otherwise        = encodeBamEntry $ trim_3' (c_trim_pred c) r'
            where r' = decodeBamEntry r

    enumInputs files >=> run                $                   -- IO ()
        joinI $ decompressBgzf              $                   -- Iteratee ByteString   IO ()
        joinI $ decodeBam                   $ \hdr ->           -- Iteratee Block        IO ()
        joinI $ mapStream do_trim           $                   -- Iteratee [BamRaw]     IO ()
        joinI $ encodeBam hdr               $                   -- Iteratee [ByteString] IO ()
        mapChunksM_ (S.hPut stdout)                             -- Iteratee ByteString   IO ()

