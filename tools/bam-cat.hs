import Bio.Bam
import System.Environment

main = do
    [inf, outf] <- getArgs
    decodeAnyBamFile inf >=> run $ \hdr ->
        mapStream (encodeBamEntry . unpackBam) =$
        writeRawBamFile outf hdr
