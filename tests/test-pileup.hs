import Bio.Adna                     ( scalarMat )
import Bio.Bam
import Bio.Bam.Pileup
import Bio.Prelude

main :: IO ()
main = do
    bams <- getArgs
    mergeInputs combineCoordinates bams >=> run $ \_ ->
            takeWhileE (isValidRefseq . b_rname . unpackBam)            =$
            mapMaybeStream (decompose (repeat (scalarMat 1)))           =$
            pileup                                                      =$
            -- mapStreamM pick                                          =$
            skipToEof

