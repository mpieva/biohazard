{-# LANGUAGE OverloadedStrings, FlexibleContexts, BangPatterns #-}
-- | Trimming of reads as found in BAM files.  Implements trimming low
-- quality sequence from the 3' end.

module Bio.Bam.Trim  {- (trim_3', trim_3, trim_low_quality) -}  where

import Bio.Bam.Header
import Bio.Bam.Rec
import Bio.Base

import Data.Bits
import Data.List ( inits, sortOn )

import qualified Data.Vector.Fusion.Stream      as S
import qualified Data.Vector.Hybrid.Internal    as Hybrid
import qualified Data.Vector.Generic            as V
import qualified Data.Vector.Unboxed            as U

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
                            zip (inits . reverse . V.toList $ b_seq b)
                                (inits . reverse . V.toList $ b_qual b)
               in trim_3 l b

    trim_rev = let l = subtract 1 . fromIntegral . length . takeWhile (uncurry p) $
                            zip (inits . V.toList $ b_seq  b)
                                (inits . V.toList $ b_qual b)
               in trim_3 l b

trim_3 :: Int -> BamRec -> BamRec
trim_3 l b | b_flag b `testBit` 4 = trim_rev
           | otherwise            = trim_fwd
  where
    trim_fwd = let (_, cigar') = trim_back_cigar (b_cigar b) l
               in b { b_seq   = V.take (V.length (b_seq  b) - l) (b_seq  b)
                    , b_qual  = V.take (V.length (b_qual b) - l) (b_qual b)
                    , b_cigar = cigar'
                    , b_exts  = deleteE "MD" (b_exts b) }

    trim_rev = let (off, cigar') = trim_fwd_cigar (b_cigar b) l
               in b { b_seq   = V.drop l (b_seq  b)
                    , b_qual  = V.drop l (b_qual b)
                    , b_cigar = cigar'
                    , b_exts  = deleteE "MD" (b_exts b)
                    , b_pos   = b_pos b + off
                    }

trim_back_cigar, trim_fwd_cigar :: V.Vector v Cigar => v Cigar -> Int -> ( Int, v Cigar )
trim_back_cigar c l = (o, V.fromList $ reverse c') where (o,c') = sanitize_cigar . trim_cigar l $ reverse $ V.toList c
trim_fwd_cigar  c l = (o, V.fromList           c') where (o,c') = sanitize_cigar $ trim_cigar l $ V.toList c

sanitize_cigar :: (Int, [Cigar]) -> (Int, [Cigar])
sanitize_cigar (o, [        ])                          = (o, [])
sanitize_cigar (o, (op:*l):xs) | op == Pad              = sanitize_cigar (o,xs)         -- del P
                               | op == Del || op == Nop = sanitize_cigar (o + l, xs)    -- adjust D,N
                               | op == Ins              = (o, (SMa :* l):xs)            -- I --> S
                               | otherwise              = (o, (op :* l):xs)             -- rest is fine

trim_cigar :: Int -> [Cigar] -> (Int, [Cigar])
trim_cigar 0 cs = (0, cs)
trim_cigar _ [] = (0, [])
trim_cigar l ((op:*ll):cs) | bad_op op = let (o,cs') = trim_cigar l cs in (o + reflen op ll, cs')
                           | otherwise = case l `compare` ll of
    LT -> (reflen op  l, (op :* (ll-l)):cs)
    EQ -> (reflen op ll,                cs)
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


-- | Overlap-merging of read pairs.  We shall compute the likelihood
-- for every possible overlap, then select the most likely one (unless it
-- looks completely random), compute a quality from the second best
-- merge, then merge and clamp the quality accordingly.
-- (We could try looking for chimaera after completing the merge, if
-- only we knew which ones to expect?)
--
-- Two reads go in, with two adapter lists.  We return 'Nothing' if all
-- merges looked mostly random.  Else we return the two original reads,
-- flagged as 'eflagVestigial' *and* the merged version, flagged as
-- 'eflagMerged' and optionally 'eflagTrimmed'.  All reads contain the
-- computed qualities (in YM and YN), which we also return.
--
-- The merging automatically limits quality scores some of the time.  We
-- additionally impose a hard limit of 63 to avoid difficulties
-- representing the result, and even that is ridiculous.  Sane people
-- would further limit the returned quality!  (In practice, map quality
-- later imposes a limit anyway, so no worries...)

merge_overlap :: BamRec -> [ U.Vector Nucleotides ]
              -> BamRec -> [ U.Vector Nucleotides ]
              -> Maybe ( BamRec, BamRec, BamRec, Int, Int )
merge_overlap r1 ads1 r2 ads2 =
    case possible_merges of
        [                             ] -> Nothing
        (score,  len) : [             ] -> result len score  plain_score
        (score1, len) : (score2, _) : _ -> result len score1 score2
  where
    rq1 = Hybrid.V (b_seq r1) (b_qual r1)
    rq2 = Hybrid.V (b_seq r2) (b_qual r2)

    -- the "merge" score if there is no overlap
    plain_score = 6 * fromIntegral (V.length rq1 + V.length rq2)

    possible_merges = sortOn fst $ [ ( merge_score ads1 ads2 rq1 rq2 l, l )
                                   | l <- [0 .. V.length rq1 + V.length rq2 - 1] ]

    flag_vestigial    br = br { b_exts = updateE "FF" (Int $ extAsInt 0 "FF" br .|. eflagVestigial) $ b_exts br }
    store_quals s1 s2 br = br { b_exts = updateE "YM" (Int $ s2          - s1) $
                                         updateE "YN" (Int $ plain_score - s1) $ b_exts br }

    result l s1 s2 = Just ( store_quals s1 s2 $ flag_vestigial r1
                          , store_quals s1 s2 $ flag_vestigial r2
                          , store_quals s1 s2 $ merged_read l (fromIntegral . min 63 $ s2-s1)
                          , s2 - s1, plain_score - s1 )

    -- How do we actually trim the read?
    -- - compute new sequence, new qualities
    -- - limit qualities to YM
    merged_read l qmax | V.length merged_seq == l
        = nullBamRec {
            b_qname = b_qname r1,
            b_flag  = flagUnmapped .|. complement pair_flags .&. b_flag r1,
            b_seq   = V.unstream $ S.map fst $ V.stream merged_seq,
            b_qual  = V.unstream $ S.map snd $ V.stream merged_seq,
            b_exts  = let ff = if l < V.length (b_seq r1) then eflagTrimmed else 0
                      in updateE "FF" (Int $ extAsInt 0 "FF" r1 .|. eflagMerged .|. ff) $ b_exts r1 }
      where
        merged_seq = V.concat
                [ V.take (l - V.length rq2) rq1
                , merge_seqs qmax        (V.take l $ V.drop (l - V.length rq2) rq1)
                             (V.reverse $ V.take l $ V.drop (l - V.length rq1) rq2)
                , V.reverse $ V.take (l - V.length rq1) rq2 ]

    pair_flags = flagPaired.|.flagProperlyPaired.|.flagMateUnmapped.|.flagMateReversed.|.flagFirstMate.|.flagSecondMate

    merge_seqs qmax = V.zipWith $ \(!n1,!(Q q1)) (!n2,!(Q q2)) ->
            if     n1 == n2 then (n1, Q $ min qmax (q1 + q2))
            else if q1 > q2 then (n1, Q $           q1 - q2 )
            else                 (n2, Q $           q2 - q1 )


-- | Trimming for a single read:  we need one adapter only (the one coming
-- /after/ the read), here provided as a list of options, and then we
-- merge with an empty second read.  Results in up to two reads (the
-- original, possibly flagged, and the trimmed one, definitely flagged,
-- and two qualities).
-- trim_adapter :: BamRec -> [???] -> Maybe (BamRec, BamRec, Qual, Qual)



-- | For merging, we don't need the complete adapters (length around 70!),
-- only a sufficient prefix.  Taking only the more-or-less constant
-- part (length around 30), there aren't all that many different
-- adapters in the world.  To deal with pretty much every library, we
-- only need the following forward adapters, which will be the default
-- (defined here in the direction they would be sequenced in):  Genomic
-- R2, Multiplex R2, Fraft P7.

default_fwd_adapters :: [ U.Vector Nucleotides ]
default_fwd_adapters = map (U.fromList. map toNucleotides)
         [ {- Genomic R2   -}  "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG"
         , {- Multiplex R2 -}  "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
         , {- Graft P7     -}  "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG" ]

-- | Like 'default_rev_adapters', these are the few adapters needed for
-- the reverse read (defined in the direction they would be sequenced in
-- as part of the second read):  Genomic R1, CL 72.

default_rev_adapters :: [ U.Vector Nucleotides ]
default_rev_adapters = map (U.fromList. map toNucleotides)
         [ {- Genomic_R1   -}  "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
         , {- CL72         -}  "GGAAGAGCGTCGTGTAGGGAAAGAGTGT" ]

-- We need to compute the likelihood of a read pair given an assumed
-- insert length.  The likelihood of the first read is the likelihood of
-- a match with the adapter where it overlaps the 3' adapter, elsewhere
-- it's 1/4 per position.  The likelihood of the second read is the
-- likelihood of a match with the adapter where it overlaps the adapter,
-- the likehood of a read-read match where it overlaps read one, 1/4 per
-- position elsewhere.  (Yes, this ignores base composition.  It doesn't
-- matter enough.)

merge_score
    :: ( V.Vector v Nucleotides, V.Vector u (Nucleotides, Qual) )
    => [ v Nucleotides ]         -- 3' adapters as they appear in the first read
    -> [ v Nucleotides ]         -- 5' adapters as they appear in the second read
    -> u (Nucleotides, Qual)            -- first read
    -> u (Nucleotides, Qual)           -- second read
    -> Int                          -- assumed insert length
    -> Int                          -- score (roughly sum of qualities at mismatches)
merge_score fwd_adapters rev_adapters read1 read2 l
    =   6 * fromIntegral (l `min` V.length read1)                                           -- read1, part before adapter
      + 6 * fromIntegral (max 0 (l - V.length read1))                                       -- read2, part before overlap

      + minimum [ match_adapter (V.drop l read1) fwd_ad                                     -- read1, match with forward adapter
                + 6 * fromIntegral (max 0 (V.length read1 - V.length fwd_ad - l))           -- read1, part after (known) adapter
                | fwd_ad <- fwd_adapters ]

      + minimum [ match_adapter (V.drop l read2) rev_ad                                     -- read2, match with reverse adapter
                + 6 * fromIntegral (max 0 (V.length read2 - V.length rev_ad - l))           -- read2, part after (known) adapter
                | rev_ad <- rev_adapters ]

      + match_reads (V.take l $ V.drop (l - V.length read2) read1)
                    (V.take l $ V.drop (l - V.length read1) read2)                          -- read2, overlap with read1
  where
    -- match_adapter :: u (Nucleotides, Qual) -> v Nucleotides -> Double
    match_adapter rd ad = S.foldl' (+) 0 $
                          S.zipWith (\(!n, Q !q) m -> if n == m then 0 else min 25 (fromIntegral q))
                                    (V.stream rd) (V.stream ad)

    -- match_reads :: u (Nucleotides, Qual) -> u (Nucleotides, Qual) -> Double
    match_reads rd1 rd2 = S.foldl' (+) 0 $
                          S.zipWith (\(!n1, Q !q1) (!n2, Q !q2) -> if n1 == compls n2 then 0 else fromIntegral $ min q1 q2)
                                    (V.stream rd1) (V.streamR rd2)
