import Bio.Bam
import System.Environment
import System.IO
import Data.ByteString.Char8 ( hPut )

main = do
    [m, inf, outf] <- getArgs
    case m of
        "P" -> decodeAnyBamFile inf >=> run $ \hdr ->
                    mapStream (encodeBamEntry . unpackBam) =$
                    writeRawBamFile outf hdr

        "R" -> decodeAnyBamFile inf >=> run $ \hdr ->
                    writeRawBamFile outf hdr

        "B" -> withFile outf WriteMode $ \hdl ->
                    decodeAnyBamFile inf >=> run $ \hdr ->
                        mapStream unpackBam =$
                        encodeBamWith2 6 hdr =$
                        mapChunksM_ (hPut hdl)

        "r" -> withFile outf WriteMode $ \hdl ->
                    decodeAnyBamFile inf >=> run $ \hdr ->
                        encodeBamRawWith2 6 hdr =$
                        mapChunksM_ (hPut hdl)
