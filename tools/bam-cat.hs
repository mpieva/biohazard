-- Silly benchmarking program.  Reads bam and writes bam, maybe in
-- different ways.

import Bio.Bam
import System.Environment

main = do
    [fi,fo] <- getArgs
    decodeAnyBamFile fi >=> run $ \hdr ->
        mapStream unpackBam {-decodeBamEntry-} =$
        takeStream 100000 =$
        mapStream encodeBamEntry =$
        writeRawBamFile fo hdr

    
