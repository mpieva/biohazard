{-# LANGUAGE OverloadedStrings #-}
-- | This module contains stuff relating to conventions local to MPI
-- EVAN.  The code is needed regularly, but it can be harmful when
-- applied to BAM files that follow different conventions.  Most
-- importantly, no program should call these functions by default.

module Bio.Bam.Evan where

import Bio.Bam.Header
import Bio.Bam.Rec
import Data.Bits
import Prelude

import qualified Data.ByteString.Char8 as S

-- | Fixes abuse of flags valued 0x800 and 0x1000.  We used them for
-- low quality and low complexity, but they have since been redefined.
-- If set, we clear them and store them into the ZD field.  Also fixes
-- abuse of the combination of the paired, 1st mate and 2nd mate flags
-- used to indicate merging or trimming.  These are canonicalized and
-- stored into the FF field.  This function is unsafe on BAM files of
-- unclear origin!
fixupFlagAbuse :: BamRec -> BamRec
fixupFlagAbuse b =
    (if b_flag b .&. flag_low_quality /= 0 then setQualFlag 'Q' else id) $          -- low qual, new convention
    (if b_flag b .&. flag_low_complexity /= 0 then setQualFlag 'C' else id) $       -- low complexity, new convention
    b { b_flag = cleaned_flags, b_exts = cleaned_exts }
  where
        -- removes old flag abuse
        flags' = b_flag b .&. complement (flag_low_quality .|. flag_low_complexity)
        cleaned_flags | flags' .&. flagPaired == 0 = flags' .&. complement (flagFirstMate .|. flagSecondMate)
                      | otherwise                  = flags'

        flag_low_quality    =  0x800
        flag_low_complexity = 0x1000

        -- merged & trimmed from old flag abuse
        is_merged  = flags' .&. (flagPaired .|. flagFirstMate .|. flagSecondMate) == flagFirstMate .|. flagSecondMate
        is_trimmed = flags' .&. (flagPaired .|. flagFirstMate .|. flagSecondMate) == flagSecondMate
        newflags = (if is_merged then eflagMerged else 0) .|. (if is_trimmed then eflagTrimmed else 0)

        -- Extended flags, renamed to avoid collision with BWA Goes like this:  if FF is there, use
        -- it.  Else check if XF is there _and_is_numeric_.  If so, use it and remove it, and set FF
        -- instead.  Else use 0 and leave it alone.  Note that this solves the collision with BWA,
        -- since BWA puts a character there, not an int.
        cleaned_exts = case (lookup "FF" (b_exts b), lookup "XF" (b_exts b)) of
                ( Just (Int i), _ ) -> updateE "FF" (Int (i .|. newflags))                (b_exts b)
                ( _, Just (Int i) ) -> updateE "FF" (Int (i .|. newflags)) $ deleteE "XF" (b_exts b)
                _ | newflags /= 0   -> updateE "FF" (Int        newflags )                (b_exts b)
                  | otherwise       ->                                                     b_exts b


-- | Fixes typical inconsistencies produced by Bwa: sometimes, 'mate unmapped' should be set, and we
-- can see it, because we match the mate's coordinates.  Sometimes 'properly paired' should not be
-- set, because one mate in unmapped.  This function is generally safe, but needs to be called only
-- on the output of affected (older?) versions of Bwa.
fixupBwaFlags :: BamRec -> BamRec
fixupBwaFlags b = b { b_flag = fixPP $ b_flag b .|. if mu then flagMateUnmapped else 0 }
  where
        -- Set "mate unmapped" if self coordinates and mate coordinates are equal, but self is
        -- paired and mapped.  (BWA forgets this flag for invalid mate alignments)
        mu = and [ isPaired b, not (isUnmapped b)
                 , isReversed b == isMateReversed b
                 , b_rname b == b_mrnm b, b_pos b == b_mpos b ]

        -- If either mate is unmapped, remove "properly paired".
        fixPP f | f .&. (flagUnmapped .|. flagMateUnmapped) == 0 = f
                | otherwise = f .&. complement flagProperlyPaired

-- | Removes syntactic warts from old read names or the read names used
-- in FastQ files.
removeWarts :: BamRec -> BamRec
removeWarts br = br { b_qname = name, b_flag = flags, b_exts = tags }
  where
    (name, flags, tags) = checkFR $ checkC $ checkSharp (b_qname br, b_flag br, b_exts br)

    checkFR (n,f,t) | "F_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagFirstMate  .|. flagPaired, t)
                    | "R_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagSecondMate .|. flagPaired, t)
                    | "M_" `S.isPrefixOf` n = checkC (S.drop 2 n, f,   insertE "FF" (Int  eflagMerged) t)
                    | "T_" `S.isPrefixOf` n = checkC (S.drop 2 n, f,   insertE "FF" (Int eflagTrimmed) t)
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


