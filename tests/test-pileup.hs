import Bio.Bam
import Bio.Bam.Pileup
import Bio.Prelude

main :: IO ()
main = do
    bams <- getArgs
    mergeInputs combineCoordinates bams >=> run                  $ \_ ->
            takeWhileE (isValidRefseq . b_rname . unpackBam)    =$
            concatMapStream (decompose $ DmgToken 0)            =$
            pileup                                              =$
            skipToEof
