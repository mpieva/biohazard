import Bio.Bam
import System.Environment
import System.IO
import Data.ByteString.Char8 ( hPut )

main :: IO ()
main = do
    [m, inf, outf] <- getArgs
    case m of
        "R" -> decodeAnyBamFile inf >=> run $ \hdr ->
                    writeBamFile outf hdr

        "B" -> withFile outf WriteMode $ \hdl ->
                    decodeAnyBamFile inf >=> run $ \hdr ->
                        mapStream unpackBam =$
                        encodeBamWith 6 hdr =$
                        mapChunksM_ (hPut hdl)

        "r" -> withFile outf WriteMode $ \hdl ->
                    decodeAnyBamFile inf >=> run $ \hdr ->
                        encodeBamWith 6 hdr =$
                        mapChunksM_ (hPut hdl)
