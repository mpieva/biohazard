{-# LANGUAGE OverloadedStrings #-}
module Bio.Bam.Evan where

-- ^ This module contains stuff relating to conventions local to MPI
-- EVAN.  The code is needed regularly, but it can be harmful when
-- applied to BAM files that follow different conventions.  Most
-- importantly, no program should call these functions by default.

import Bio.Bam.Header
import Bio.Bam.Rec
import Data.Bits

import qualified Data.ByteString.Char8 as S

-- | Fixes a BAM record for changed conventions.  This should probably
-- be split into multiple functions that fix different changed
-- conventions...   XXX
--
-- Merged and trimmed flags need different handling.  We can't put them
-- into flags without imposing a performance penalty on pretty much
-- everything.   XXX
fixup_bam_rec :: BamRec -> BamRec
fixup_bam_rec b =
    (if b_flag b .&. flagLowQuality /= 0 then setQualFlag 'Q' else id) $          -- low qual, new convention
    (if b_flag b .&. flagLowComplexity /= 0 then setQualFlag 'C' else id) $       -- low complexity, new convention
    b { b_flag = fixPP $ oflags .|. muflag .|. tflags .|. shiftL eflags 16        -- extended flags
      , b_exts = cleaned_exts }
  where
        -- removes old flag abuse
        flags' = b_flag b .&. complement (flagLowQuality .|. flagLowComplexity)
        oflags | flags' .&. flagPaired == 0 = flags' .&. complement (flagFirstMate .|. flagSecondMate)
               | otherwise                  = flags'

        -- set "mate unmapped" if self coordinates and mate coordinates are equal, but self is paired and mapped
        -- (BWA forgets this flag for invalid mate alignments)
        muflag = if mu then flagMateUnmapped else 0
        mu = and [ isPaired b, not (isUnmapped b)
                 , isReversed b == isMateReversed b
                 , b_rname b == b_mrnm b, b_pos b == b_mpos b ]

        -- merged & trimmed from old flag abuse
        is_merged  = flags' .&. (flagPaired .|. flagFirstMate .|. flagSecondMate) == flagFirstMate .|. flagSecondMate
        is_trimmed = flags' .&. (flagPaired .|. flagFirstMate .|. flagSecondMate) == flagSecondMate

        tflags = (if is_merged  then flagMerged  else 0) .|.
                 (if is_trimmed then flagTrimmed else 0)

        -- extended flags, renamed to avoid collision with BWA
        -- Goes like this:  if FF is there, use and remove it.  Else
        -- check if XF is there _and_is_numeric_.  If so, use it and
        -- remove it.  Else use 0 and leave it alone.  Note that this
        -- solves the collision with BWA, since BWA puts a character
        -- there, not an int.
        (eflags, cleaned_exts) = case (lookup "FF" (b_exts b), lookup "XF" (b_exts b)) of
                ( Just (Int i), _ ) -> (i, deleteE "FF" (b_exts b))
                ( _, Just (Int i) ) -> (i, deleteE "XF" (b_exts b))
                (       _,_       ) -> (0,               b_exts b )

        -- if either mate is unmapped, remove "properly paired"
        fixPP f | f .&. (flagUnmapped .|. flagMateUnmapped) == 0 = f
                | otherwise = f .&. complement flagProperlyPaired

        flagLowQuality    =  0x800
        flagLowComplexity = 0x1000


-- | Remove syntactic warts from old read names.
removeWarts :: BamRec -> BamRec
removeWarts br = br { b_qname = name, b_flag = flags, b_exts = tags }
  where
    (name, flags, tags) = checkFR $ checkC $ checkSharp (b_qname br, b_flag br, b_exts br)

    checkFR (n,f,t) | "F_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagFirstMate  .|. flagPaired, t)
                    | "R_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagSecondMate .|. flagPaired, t)
                    | "M_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagMerged,                    t)
                    | "T_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagTrimmed,                   t)
                    | "/1" `S.isSuffixOf` n =        ( rdrop 2 n, f .|. flagFirstMate  .|. flagPaired, t)
                    | "/2" `S.isSuffixOf` n =        ( rdrop 2 n, f .|. flagSecondMate .|. flagPaired, t)
                    | otherwise             =        (         n, f,                                   t)

    checkC (n,f,t) | "C_" `S.isPrefixOf` n  = (S.drop 2 n, f, insertE "XP" (Int (-1)) t)
                   | otherwise              = (         n, f,                         t)

    rdrop n s = S.take (S.length s - n) s

    checkSharp (n,f,t) = case S.split '#' n of [n',ts] -> (n', f, insertTags ts t)
                                               _       -> ( n, f,               t)

    insertTags ts t | S.null y  = insertE "XI" (Text ts) t
                    | otherwise = insertE "XI" (Text  x) $ insertE "XJ" (Text $ S.tail y) t
        where (x,y) = S.break (== ',') ts


