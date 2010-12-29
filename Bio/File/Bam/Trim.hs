module Bio.File.Bam.Trim ( trim_3', trim_low_quality ) where

import Bio.File.Bam

import Data.Array.Unboxed
import Data.Bits ( testBit )
import Data.Int ( Int64 )
import Data.List ( genericDrop, inits )
import Data.Word ( Word8 )
import qualified Data.ByteString.Lazy as L
import qualified Data.Map as M

-- | Trims from the 3' end of a sequence.
-- @trim_3\' p b@ trims the 3' end of the sequence in @b@ at the
-- earliest position such that @p@ evaluates to true on every suffix
-- that was trimmed off.  Note that the 3' end may be the beginning of
-- the sequence if it happens to be stored in reverse-complemented form.
-- Also note that trimming from the 3' end may not make sense for reads
-- that were constructed by merging paired end data.

trim_3' :: ([Nucleotide] -> [Word8] -> Bool) -> BamRec -> BamRec
trim_3' p b | b_flag b `testBit` 4 = trim_rev
            | otherwise            = trim_fwd
  where
    trim_fwd = let l = subtract 1 . fromIntegral . length . takeWhile (uncurry p) $
                            zip (inits . reverse . decodeSeq $ b_seq b)
                                (inits . reverse . L.unpack $ b_qual b)
                   (_, cigar') = trim_back_cigar (b_cigar b) l
               in b { b_seq = trim_back_seq (b_seq b) l
                    , b_qual = L.take (L.length (b_qual b) - l) (b_qual b)
                    , b_cigar = cigar'
                    , b_exts = M.delete 17485 (b_exts b) }

    trim_rev = let l = subtract 1 . fromIntegral . length . takeWhile (uncurry p) $
                            zip (inits . decodeSeq $ b_seq b)
                                (inits . L.unpack $ b_qual b)
                   (off, cigar') = trim_fwd_cigar (b_cigar b) l
               in b { b_seq = trim_forward_seq (b_seq b) l
                    , b_qual = L.drop l (b_qual b)
                    , b_cigar = cigar'
                    , b_exts = M.delete 17485 (b_exts b)
                    , b_pos = b_pos b + off
                    }


    
trim_back_seq :: CodedSeq -> Int64 -> CodedSeq
trim_back_seq (CodedSeq len str) l = CodedSeq (len-l) (L.take ((len-l+1) `div` 2) str)

trim_forward_seq :: CodedSeq -> Int64 -> CodedSeq
trim_forward_seq sq@(CodedSeq len str) l
    | l `mod` 2 == 0 = CodedSeq (len-l) (L.drop (l `div` 2) str)
    | otherwise      = deflateSeq . genericDrop l . inflateSeq $ sq

trim_back_cigar, trim_fwd_cigar :: CodedCigar -> Int64 -> ( Int, CodedCigar )
trim_back_cigar (CodedCigar c) l = (o, packCigar $ reverse c') 
  where (o,c') = sanitize_cigar . trim_cigar l . reverse $ elems c

trim_fwd_cigar  (CodedCigar c) l = (o, packCigar           c')
  where (o,c') = sanitize_cigar . trim_cigar l           $ elems c

sanitize_cigar :: (Int, [Int]) -> (Int, [Int])
sanitize_cigar (o, [    ])                          = (o, [])
sanitize_cigar (o, (x:xs)) | op == Pad              = sanitize_cigar (o,xs)         -- del P
                           | op == Del || op == Nop = sanitize_cigar (o + l, xs)    -- adjust D,N
                           | op == Ins              = (o, mkCigOp SMa l : xs)       -- I --> S
                           | otherwise              = (o, x:xs)                     -- rest is fine
  where op = cigOp x ; l = cigLen x

trim_cigar :: Int64 -> [Int] -> (Int, [Int])
trim_cigar 0 cs = (0, cs)
trim_cigar _ [] = (0, [])
trim_cigar l (c:cs) | bad_op c = let (o,cs') = trim_cigar l cs in (o + reflen c, cs')
                    | otherwise = case l `compare` fromIntegral (cigLen c) of
    LT -> let c' = mkCigOp (cigOp c) (cigLen c - fromIntegral l) 
          in (reflen c - reflen c', c':cs)
    EQ -> (reflen c, cs)
    GT -> let (o,cs') = trim_cigar (l - fromIntegral (cigLen c)) cs
          in (o + reflen c, cs')

  where
    reflen c' = if ref_op c' then fromIntegral $ cigLen c' else 0
    bad_op c' = let o = cigOp c' in o /= Mat && o /= Ins && o /= SMa
    ref_op c' = let o = cigOp c' in o == Mat || o == Del


-- | Trim predicate to get rid of low quality sequence.
-- @trim_low_quality q ns qs@ evaluates to true if all qualities in @qs@
-- are smaller (i.e. worse) than @q@.
trim_low_quality :: Word8 -> [Nucleotide] -> [Word8] -> Bool
trim_low_quality q = \_ qs -> all (< q) qs



