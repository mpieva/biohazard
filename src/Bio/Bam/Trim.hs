-- | Trimming of reads as found in BAM files.  Implements trimming low
-- quality sequence from the 3' end.

module Bio.Bam.Trim where

import Bio.Bam.Header
import Bio.Bam.Rec
import Bio.Bam.Rmdup        ( ECig(..), setMD, toECig )
import Bio.Iteratee
import Bio.Prelude

import Foreign.C.Types      ( CInt(..) )
import Foreign.Ptr          ( Ptr )

import qualified Data.ByteString                        as B
import qualified Data.Vector.Generic                    as V
import qualified Data.Vector.Storable                   as W

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
                   c = modMd (takeECig (V.length (b_seq  b) - l)) b
               in c { b_seq   = V.take (V.length (b_seq  c) - l) (b_seq  c)
                    , b_qual  = V.take (V.length (b_qual c) - l) (b_qual c)
                    , b_cigar = cigar'
                    , b_exts  = map (\(k,e) -> case e of
                                        Text t | k `elem` trim_set
                                          -> (k, Text (B.take (B.length t - l) t))
                                        _ -> (k,e)
                                    ) (b_exts c) }

    trim_rev = let (off, cigar') = trim_fwd_cigar (b_cigar b) l
                   c = modMd (dropECig l) b
               in c { b_seq   = V.drop l (b_seq  c)
                    , b_qual  = V.drop l (b_qual c)
                    , b_pos   = b_pos c + off
                    , b_cigar = cigar'
                    , b_exts  = map (\(k,e) -> case e of
                                        Text t | k `elem` trim_set
                                          -> (k, Text (B.drop l t))
                                        _ -> (k,e)
                                    ) (b_exts c) }

    trim_set = ["BQ","CQ","CS","E2","OQ","U2"]

    modMd :: (ECig -> ECig) -> BamRec -> BamRec
    modMd f br = maybe br (setMD br . f . toECig (b_cigar br)) (getMd br)

    endOf :: ECig -> ECig
    endOf  WithMD     = WithMD
    endOf  WithoutMD  = WithoutMD
    endOf (Mat' _ es) = endOf es
    endOf (Ins' _ es) = endOf es
    endOf (SMa' _ es) = endOf es
    endOf (Rep' _ es) = endOf es
    endOf (Del' _ es) = endOf es
    endOf (Nop' _ es) = endOf es
    endOf (HMa' _ es) = endOf es
    endOf (Pad' _ es) = endOf es

    takeECig :: Int -> ECig -> ECig
    takeECig 0  es          = endOf es
    takeECig _  WithMD      = WithMD
    takeECig _  WithoutMD   = WithoutMD
    takeECig n (Mat' m  es) = Mat' n  $ if n > m then takeECig (n-m) es else WithMD
    takeECig n (Ins' m  es) = Ins' n  $ if n > m then takeECig (n-m) es else WithMD
    takeECig n (SMa' m  es) = SMa' n  $ if n > m then takeECig (n-m) es else WithMD
    takeECig n (Rep' ns es) = Rep' ns $ takeECig (n-1) es
    takeECig n (Del' ns es) = Del' ns $ takeECig n es
    takeECig n (Nop' m  es) = Nop' m  $ takeECig n es
    takeECig n (HMa' m  es) = HMa' m  $ takeECig n es
    takeECig n (Pad' m  es) = Pad' m  $ takeECig n es

    dropECig :: Int -> ECig -> ECig
    dropECig 0  es         = es
    dropECig _  WithMD     = WithMD
    dropECig _  WithoutMD  = WithoutMD
    dropECig n (Mat' m es) = if n > m then dropECig (n-m) es else Mat' n WithMD
    dropECig n (Ins' m es) = if n > m then dropECig (n-m) es else Ins' n WithMD
    dropECig n (SMa' m es) = if n > m then dropECig (n-m) es else SMa' n WithMD
    dropECig n (Rep' _ es) = dropECig (n-1) es
    dropECig n (Del' _ es) = dropECig n es
    dropECig n (Nop' _ es) = dropECig n es
    dropECig n (HMa' _ es) = dropECig n es
    dropECig n (Pad' _ es) = dropECig n es


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

merge_overlap :: BamRec -> [ W.Vector Nucleotides ]
              -> BamRec -> [ W.Vector Nucleotides ]
              -> Maybe ( BamRec, BamRec, BamRec, Int, Int )
merge_overlap r1 ads1 r2 ads2
    | V.null (b_seq r1) && V.null (b_seq r2) = Nothing
    | otherwise                              = result mlen score1 score2
  where
    -- the "merge" score if there is no overlap
    plain_score = 6 * fromIntegral (len_r1 + len_r2)

    len_r1    = V.length  $ b_seq  r1
    len_r2    = V.length  $ b_seq  r2

    b_seq_r1  = V.convert $ b_seq  r1
    b_seq_r2  = V.convert $ b_seq  r2
    b_qual_r1 = V.convert $ b_qual r1
    b_qual_r2 = V.convert $ b_qual r2

    (score1, mlen, score2) = twoMins plain_score (len_r1 + len_r2) $
                             merge_score ads1 ads2 b_seq_r1 b_qual_r1 b_seq_r2 b_qual_r2

    flag_vestigial    br = br { b_exts = updateE "FF" (Int $ extAsInt 0 "FF" br .|. eflagVestigial) $ b_exts br }
    store_quals s1 s2 br = br { b_exts = updateE "YM" (Int $ s2          - s1) $
                                         updateE "YN" (Int $ plain_score - s1) $ b_exts br }

    result l s1 s2 = Just ( store_quals s1 s2 $ flag_vestigial r1
                          , store_quals s1 s2 $ flag_vestigial r2
                          , store_quals s1 s2 $ merged_read l (fromIntegral . min 63 $ s2-s1)
                          , s2 - s1, plain_score - s1 )

    merged_read l qmax
        | V.length merged_seq /= l = error $ "Logic error in merged_read: " ++ show (V.length merged_seq, l)
        | otherwise = nullBamRec {
                b_qname = b_qname r1,
                b_flag  = flagUnmapped .|. complement pair_flags .&. b_flag r1,
                b_seq   = merged_seq,
                b_qual  = merged_qual,
                b_exts  = let ff = if l < len_r1 then eflagTrimmed else 0
                          in updateE "FF" (Int $ extAsInt 0 "FF" r1 .|. eflagMerged .|. ff) $ b_exts r1 }
      where
        merged_seq = V.convert $ V.concat
                [ V.take (l - len_r2) (b_seq_r1)
                , merge_seqs             (V.take l $ V.drop (l - len_r2) b_seq_r1)
                                         (V.take l $ V.drop (l - len_r2) b_qual_r1)
                             (V.reverse $ V.take l $ V.drop (l - len_r1) b_seq_r2)
                             (V.reverse $ V.take l $ V.drop (l - len_r1) b_qual_r2)
                , V.reverse $ V.take (l - len_r1) b_seq_r2 ]

        merged_qual = V.convert $ V.concat
                [ V.take (l - len_r2) (b_qual_r1)
                , merge_quals qmax        (V.take l $ V.drop (l - len_r2) b_seq_r1)
                                          (V.take l $ V.drop (l - len_r2) b_qual_r1)
                              (V.reverse $ V.take l $ V.drop (l - len_r1) b_seq_r2)
                              (V.reverse $ V.take l $ V.drop (l - len_r1) b_qual_r2)
                , V.reverse $ V.take (l - len_r1) b_qual_r2 ]

    pair_flags = flagPaired.|.flagProperlyPaired.|.flagMateUnmapped.|.flagMateReversed.|.flagFirstMate.|.flagSecondMate

    merge_seqs v1 v2 v3 v4 = V.zipWith4 zz v1 v2 v3 v4
      where
        zz !n1 (Q !q1) !n2 (Q !q2) | n1 == compls n2 =        n1
                                   | q1 > q2         =        n1
                                   | otherwise       = compls n2

    merge_quals qmax v1 v2 v3 v4 = V.zipWith4 zz v1 v2 v3 v4
      where
        zz !n1 (Q !q1) !n2 (Q !q2) | n1 == compls n2 = Q $ min qmax (q1 + q2)
                                   | q1 > q2         = Q $           q1 - q2
                                   | otherwise       = Q $           q2 - q1

-- | Trimming for a single read:  we need one adapter only (the one coming
-- /after/ the read), here provided as a list of options, and then we
-- merge with an empty second read.  Results in up to two reads (the
-- original, possibly flagged, and the trimmed one, definitely flagged,
-- and two qualities).
trim_adapter :: BamRec -> [ W.Vector Nucleotides ] -> Maybe ( BamRec, BamRec, Int, Int )
trim_adapter r1 ads1
    | V.null (b_seq r1) = Nothing
    | otherwise         = result mlen score1 score2
  where
    -- the "merge" score if there is no trimming
    plain_score = 6 * fromIntegral (V.length (b_seq r1))

    b_seq_r1 = V.convert $ b_seq r1
    b_qual_r1 = V.convert $ b_qual r1

    (score1, mlen, score2) = twoMins plain_score (V.length (b_seq r1)) $
                             merge_score ads1 [V.empty] b_seq_r1 b_qual_r1 V.empty V.empty

    flag_vestigial    br = br { b_exts = updateE "FF" (Int $ extAsInt 0 "FF" br .|. eflagVestigial) $ b_exts br }
    store_quals s1 s2 br = br { b_exts = updateE "YM" (Int $ s2          - s1) $
                                         updateE "YN" (Int $ plain_score - s1) $ b_exts br }

    result l s1 s2 = Just ( store_quals s1 s2 $ flag_vestigial r1
                          , store_quals s1 s2 $ trimmed_read l
                          , s2 - s1, plain_score - s1 )

    trimmed_read l = nullBamRec {
            b_qname = b_qname r1,
            b_flag  = flagUnmapped .|. b_flag r1,
            b_seq   = V.take l $ b_seq  r1,
            b_qual  = V.take l $ b_qual r1,
            b_exts  = updateE "FF" (Int $ extAsInt 0 "FF" r1 .|. eflagTrimmed) $ b_exts r1 }


-- | For merging, we don't need the complete adapters (length around 70!),
-- only a sufficient prefix.  Taking only the more-or-less constant
-- part (length around 30), there aren't all that many different
-- adapters in the world.  To deal with pretty much every library, we
-- only need the following forward adapters, which will be the default
-- (defined here in the direction they would be sequenced in):  Genomic
-- R2, Multiplex R2, Fraft P7.

default_fwd_adapters :: [ W.Vector Nucleotides ]
default_fwd_adapters = map (W.fromList. map toNucleotides)
         [ {- Genomic R2   -}  "AGATCGGAAGAGCGGTTCAG"
         , {- Multiplex R2 -}  "AGATCGGAAGAGCACACGTC"
         , {- Graft P7     -}  "AGATCGGAAGAGCTCGTATG" ]

-- | Like 'default_rev_adapters', these are the few adapters needed for
-- the reverse read (defined in the direction they would be sequenced in
-- as part of the second read):  Genomic R1, CL 72.

default_rev_adapters :: [ W.Vector Nucleotides ]
default_rev_adapters = map (W.fromList. map toNucleotides)
         [ {- Genomic_R1   -}  "AGATCGGAAGAGCGTCGTGT"
         , {- CL72         -}  "GGAAGAGCGTCGTGTAGGGA" ]

-- We need to compute the likelihood of a read pair given an assumed
-- insert length.  The likelihood of the first read is the likelihood of
-- a match with the adapter where it overlaps the 3' adapter, elsewhere
-- it's 1/4 per position.  The likelihood of the second read is the
-- likelihood of a match with the adapter where it overlaps the adapter,
-- the likehood of a read-read match where it overlaps read one, 1/4 per
-- position elsewhere.  (Yes, this ignores base composition.  It doesn't
-- matter enough.)

merge_score
    :: [ W.Vector Nucleotides ]                 -- 3' adapters as they appear in the first read
    -> [ W.Vector Nucleotides ]                 -- 5' adapters as they appear in the second read
    -> W.Vector Nucleotides -> W.Vector Qual    -- first read, qual
    -> W.Vector Nucleotides -> W.Vector Qual    -- second read, qual
    -> Int                                      -- assumed insert length
    -> Int                                      -- score (roughly sum of qualities at mismatches)
merge_score fwd_adapters rev_adapters !read1 !qual1 !read2 !qual2 !l
    =   6 * fromIntegral (l `min` V.length read1)                                           -- read1, part before adapter
      + 6 * fromIntegral (max 0 (l - V.length read1))                                       -- read2, part before overlap

      + foldl' (\acc fwd_ad -> min acc
                    (match_adapter l read1 qual1 fwd_ad +                                   -- read1, match with forward adapter
                     6 * fromIntegral (max 0 (V.length read1 - V.length fwd_ad - l)))       -- read1, part after (known) adapter
               ) maxBound fwd_adapters

      + foldl' (\acc rev_ad -> min acc
                    (match_adapter l read2 qual2 rev_ad +                                   -- read2, match with reverse adapter
                     6 * fromIntegral (max 0 (V.length read2 - V.length rev_ad - l)))       -- read2, part after (known) adapter
               ) maxBound rev_adapters

      + match_reads l read1 qual1 read2 qual2

{-# INLINE match_adapter #-}
match_adapter :: Int -> W.Vector Nucleotides -> W.Vector Qual -> W.Vector Nucleotides -> Int
match_adapter !off !rd !qs !ad
    | V.length rd /= V.length qs = error "read/qual length mismatch"
    | efflength <= 0             = 0
    | otherwise
        = fromIntegral . unsafePerformIO $
          W.unsafeWith rd $ \p_rd ->
          W.unsafeWith qs $ \p_qs ->
          W.unsafeWith ad $ \p_ad ->
          prim_match_ad (fromIntegral off)
                        (fromIntegral efflength)
                        p_rd p_qs p_ad
  where
    !efflength =  (V.length rd - off) `min` V.length ad

foreign import ccall unsafe "prim_match_ad"
    prim_match_ad :: CInt -> CInt
                  -> Ptr Nucleotides -> Ptr Qual
                  -> Ptr Nucleotides -> IO CInt


-- | Computes overlap score for two reads (with qualities) assuming an
-- insert length.
{-# INLINE match_reads #-}
match_reads :: Int -> W.Vector Nucleotides -> W.Vector Qual -> W.Vector Nucleotides -> W.Vector Qual -> Int
match_reads !l !rd1 !qs1 !rd2 !qs2
    | V.length rd1 /= V.length qs1 || V.length rd2 /= V.length qs2 = error "read/qual length mismatch"
    | efflength <= 0                                               = 0
    | otherwise
        = fromIntegral . unsafePerformIO $
          W.unsafeWith rd1 $ \p_rd1 ->
          W.unsafeWith qs1 $ \p_qs1 ->
          W.unsafeWith rd2 $ \p_rd2 ->
          W.unsafeWith qs2 $ \p_qs2 ->
          prim_match_reads (fromIntegral minidx1)
                           (fromIntegral maxidx2)
                           (fromIntegral efflength)
                           p_rd1 p_qs1 p_rd2 p_qs2
  where
    -- vec1, forward
    !minidx1 = (l - V.length rd2) `max` 0
    -- vec2, backward
    !maxidx2 = l `min` V.length rd2
    -- effective length
    !efflength = ((V.length rd1 + V.length rd2 - l) `min` l) `max` 0


foreign import ccall unsafe "prim_match_reads"
    prim_match_reads :: CInt -> CInt -> CInt
                     -> Ptr Nucleotides -> Ptr Qual
                     -> Ptr Nucleotides -> Ptr Qual -> IO CInt


{-# INLINE twoMins #-}
twoMins :: (Bounded a, Ord a) => a -> Int -> (Int -> a) -> (a,Int,a)
twoMins a0 imax f = go a0 0 maxBound 0 0
  where
    go !m1 !i1 !m2 !i2 !i
        | i == imax = (m1,i1,m2)
        | otherwise =
            case f i of
                x | x < m1    -> go  x  i m1 i1 (i+1)
                  | x < m2    -> go m1 i1  x  i (i+1)
                  | otherwise -> go m1 i1 m2 i2 (i+1)


mergeTrimBam :: Monad m => [W.Vector Nucleotides] -> [W.Vector Nucleotides] -> Enumeratee [BamRec] [BamRec] m a
mergeTrimBam fwd_ads rev_ads = convStream go
  where
    go = do r1 <- headStream
            if isPaired r1
              then tryHead >>= go2 r1
              else case trim_adapter r1 fwd_ads of
                    Nothing                -> return [r1]
                    Just (r1',r1t,_q1,_q2) -> return [r1t,r1']

    go2 r1  Nothing  = error $ "Lone mate found: " ++ show (b_qname r1)
    go2 r1 (Just r2) = case merge_overlap r1 fwd_ads r2 rev_ads of
                    Nothing                   -> return [r1,r2]
                    Just (r1',r2',rm,_q1,_q2) -> return [rm,r1',r2']

