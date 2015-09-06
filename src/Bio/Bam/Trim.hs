-- | Trimming of reads as found in BAM files.  Implements trimming low
-- quality sequence from the 3' end.

module Bio.Bam.Trim ( trim_3', trim_3, trim_low_quality ) where

import Bio.Bam.Header
import Bio.Bam.Rec
import Bio.Base

import Data.Bits ( testBit )
import Data.List ( inits )
import qualified Data.ByteString     as S
import qualified Data.Map            as M
import qualified Data.Vector.Generic as V

-- | Trims from the 3' end of a sequence.
-- @trim_3\' p b@ trims the 3' end of the sequence in @b@ at the
-- earliest position such that @p@ evaluates to true on every suffix
-- that was trimmed off.  Note that the 3' end may be the beginning of
-- the sequence if it happens to be stored in reverse-complemented form.
-- Also note that trimming from the 3' end may not make sense for reads
-- that were constructed by merging paired end data (but we cannot take
-- care of that here).  Further note that trimming may break dependent
-- information, notably the "mate" information of the mate and many
-- optional fields.
--
-- TODO: The MD field is currently removed.  It should be repaired
-- instead.  Many other fields should be trimmed if present.

trim_3' :: ([Nucleotides] -> [Qual] -> Bool) -> BamRec -> BamRec
trim_3' p b | b_flag b `testBit` 4 = trim_rev
            | otherwise            = trim_fwd
  where
    trim_fwd = let l = subtract 1 . fromIntegral . length . takeWhile (uncurry p) $
                            zip (inits . reverse .         V.toList $ b_seq b)
                                (inits . reverse . map Q . S.unpack $ b_qual b)
               in trim_3 l b

    trim_rev = let l = subtract 1 . fromIntegral . length . takeWhile (uncurry p) $
                            zip (inits .         V.toList $ b_seq  b)
                                (inits . map Q . S.unpack $ b_qual b)
               in trim_3 l b

trim_3 :: Int -> BamRec -> BamRec
trim_3 l b | b_flag b `testBit` 4 = trim_rev
           | otherwise            = trim_fwd
  where
    trim_fwd = let (_, cigar') = trim_back_cigar (b_cigar b) l
               in b { b_seq   = V.take (V.length (b_seq  b) - l) (b_seq  b)
                    , b_qual  = S.take (S.length (b_qual b) - l) (b_qual b)
                    , b_cigar = cigar'
                    , b_exts  = filter ((/=) "MD" . fst) (b_exts b) }

    trim_rev = let (off, cigar') = trim_fwd_cigar (b_cigar b) l
               in b { b_seq   = V.drop l (b_seq  b)
                    , b_qual  = S.drop l (b_qual b)
                    , b_cigar = cigar'
                    , b_exts  = filter ((/=) "MD" . fst) (b_exts b)
                    , b_pos   = b_pos b + off
                    }

trim_back_cigar, trim_fwd_cigar :: Cigar -> Int -> ( Int, Cigar )
trim_back_cigar (Cigar c) l = (o, Cigar $ reverse c') where (o,c') = sanitize_cigar . trim_cigar l $ reverse c
trim_fwd_cigar  (Cigar c) l = (o, Cigar           c') where (o,c') = sanitize_cigar $ trim_cigar l c

sanitize_cigar :: (Int, [(CigOp,Int)]) -> (Int, [(CigOp,Int)])
sanitize_cigar (o, [       ])                          = (o, [])
sanitize_cigar (o, (op,l):xs) | op == Pad              = sanitize_cigar (o,xs)         -- del P
                              | op == Del || op == Nop = sanitize_cigar (o + l, xs)    -- adjust D,N
                              | op == Ins              = (o, (SMa,l) : xs)             -- I --> S
                              | otherwise              = (o, (op,l):xs)                -- rest is fine

trim_cigar :: Int -> [(CigOp,Int)] -> (Int, [(CigOp,Int)])
trim_cigar 0 cs = (0, cs)
trim_cigar _ [] = (0, [])
trim_cigar l ((op,ll):cs) | bad_op op = let (o,cs') = trim_cigar l cs in (o + reflen op ll, cs')
                          | otherwise = case l `compare` ll of
    LT -> (reflen op  l, (op,ll-l):cs)
    EQ -> (reflen op ll,           cs)
    GT -> let (o,cs') = trim_cigar (l - ll) cs in (o + reflen op ll, cs')

  where
    reflen op' = if ref_op op' then id else const 0
    bad_op o = o /= Mat && o /= Ins && o /= SMa
    ref_op o = o == Mat || o == Del


-- | Trim predicate to get rid of low quality sequence.
-- @trim_low_quality q ns qs@ evaluates to true if all qualities in @qs@
-- are smaller (i.e. worse) than @q@.
trim_low_quality :: Qual -> a -> [Qual] -> Bool
trim_low_quality q = const $ all (< q)



