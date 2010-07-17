module Bio.File.Bam.Mpi (
    module Bio.File.Bam,

    bamFlagQualityFiltered,
    bamFlagComplexityFiltered,
    bamFilterFlags,
    isFirstMate,
    isSecondMate,
    isAdapterTrimmed
) where

import Bio.File.Bam
import Data.Bits

-- * Bam functions depending on MPI conventions

bamFlagQualityFiltered, bamFlagComplexityFiltered :: Int
bamFlagQualityFiltered = 0x800
bamFlagComplexityFiltered = 0x1000

-- | get all the filter flags
-- Filter flags are the standard "fails QC", "is duplicate" and our
-- extensions "fails quality" and "fails complexity".
bamFilterFlags :: BamRec -> Int
bamFilterFlags = (.&. mask) . b_flag  
  where mask = bamFlagQualityFiltered .|. bamFlagComplexityFiltered .|. bamFlagFailsQC

-- | tests if a record is a "first mate"
-- Returns true if the read is flagged as "paired" and "first mate".
-- Does not return true for single reads.
isFirstMate :: BamRec -> Bool
isFirstMate = (== good) . (.&. mask) . b_flag
  where good = bamFlagPaired .|. bamFlagFirstMate
        mask = good .|. bamFlagSecondMate

-- | tests if a record is a "second mate"
-- Returns true if the read is flagged as "paired" and "second mate".
-- Does not return true for single reads, even if they are adapter
-- trimmed.
isSecondMate :: BamRec -> Bool
isSecondMate = (== good) . (.&. mask) . b_flag
  where good = bamFlagPaired .|. bamFlagSecondMate
        mask = good .|. bamFlagFirstMate

-- | tests if a read is adapter trimmed
-- Returns true for adapter trimmed single reads, does not get confused
-- by paired reads.
isAdapterTrimmed :: BamRec -> Bool
isAdapterTrimmed = (== good) . (.&. mask) . b_flag
  where good = bamFlagFirstMate .|. bamFlagSecondMate
        mask = bamFlagPaired .|. good

