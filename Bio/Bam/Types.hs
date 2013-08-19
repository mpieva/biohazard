{-# LANGUAGE OverloadedStrings #-}

module Bio.Bam.Types (
        Cigar(..),
        CigOp(..),
        cigarToAlnLen
    ) where

import Data.Ix (Ix)
import Data.ByteString.Char8 (index)

-- ^ Types for representing the BAM data model that are useful in more
-- than one place.

-- | Cigar line in BAM coding
-- Bam encodes an operation and a length into a single integer, we keep
-- those integers in an array.
newtype Cigar = Cigar { unCigar :: [(CigOp, Int)] }

data CigOp = Mat | Ins | Del | Nop | SMa | HMa | Pad
    deriving ( Eq, Ord, Enum, Show, Bounded, Ix )

instance Show Cigar where
    show (Cigar cs) = concat [ shows l (toChr op) | (op,l) <- cs ]
      where toChr = (:[]) . index "MIDNSHP=X" . fromEnum


-- | extracts the aligned length from a cigar line
-- This gives the length of an alignment as measured on the reference,
-- which is different from the length on the query or the length of the
-- alignment.
cigarToAlnLen :: Cigar -> Int
cigarToAlnLen (Cigar cig) = sum $ map l cig
  where l (op,n) = if op == Mat || op == Del || op == Nop then n else 0


