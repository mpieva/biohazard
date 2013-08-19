{-# LANGUAGE ExistentialQuantification, RecordWildCards, NamedFieldPuns #-}
module Bio.File.Bam.Rmdup(
            rmdup, Collapse, cons_collapse, cheap_collapse,
            cons_collapse_keep, cheap_collapse_keep,
            check_sort
    ) where

import Bio.Bam
import Bio.Base
import Bio.File.Bam.Fastq               ( removeWarts )
import Bio.Iteratee
import Data.Array.Unboxed
import Data.Bits
import Data.List
import Data.Ord                         ( comparing )
import Data.Word                        ( Word8 )

import qualified Data.ByteString        as B
import qualified Data.ByteString.Char8  as T
import qualified Data.Iteratee          as I
import qualified Data.Map               as M
import qualified Data.Vector.Generic    as V

data Collapse = forall a . Collapse {
                    inject :: BamRaw -> a,                  -- convert BamRaw to internal data type
                    collapse :: [a] -> (Either a a,[a]),    -- cluster to consensus and stuff or representative and stuff
                    add_xp_of :: a -> a -> a,               -- modify XP when discarding a record
                    originals :: [a] -> [a],                -- treatment of the redundant original reads
                    make_singleton :: a -> a,               -- force to singleton
                    project :: a -> BamRaw }                -- get it back out

cons_collapse :: Word8 -> Collapse
cons_collapse maxq = Collapse (removeWarts . decodeBamEntry) (do_collapse maxq)
                              addXPOf (const []) make_singleton_cooked encodeBamEntry

cons_collapse_keep :: Word8 -> Collapse
cons_collapse_keep maxq = Collapse (removeWarts . decodeBamEntry) (do_collapse maxq)
                                   addXPOf (map flagDup) make_singleton_cooked encodeBamEntry
  where
    flagDup b = b { b_flag = b_flag b .|. flagDuplicate }

cheap_collapse :: Collapse
cheap_collapse = Collapse id do_cheap_collapse addXPOf' (const []) make_singleton_raw id

cheap_collapse_keep :: Collapse
cheap_collapse_keep = Collapse id do_cheap_collapse addXPOf' (map flagDup) make_singleton_raw id
  where
    flagDup b = mutateBamRaw b $ setFlag (br_flag b .|. flagDuplicate)

make_singleton_raw :: BamRaw -> BamRaw
make_singleton_raw b = mutateBamRaw b $ setFlag (br_flag b .&. complement pflags)

make_singleton_cooked :: BamRec -> BamRec
make_singleton_cooked b = b { b_flag = b_flag b .&. complement pflags }

pflags :: Int
pflags = flagPaired .|. flagProperlyPaired .|. flagMateUnmapped .|. flagMateReversed .|. flagFirstMate .|. flagSecondMate

-- | Removes duplicates from an aligned, sorted BAM stream.
--
-- The incoming stream must be sorted by coordinate, and we check for
-- violations of that assumption.  We cannot assume that length was
-- taken into account when sorting (samtools doesn't do so), so
-- duplicates may be separated by reads that start at the same position
-- but have different length or different strand.
--
-- We are looking at three different kinds of reads:  paired reads, true
-- single ended reads, merged or trimmed reads.  They are somewhat
-- different, but here's the situation if we wanted to treat them
-- separately.  These conditions define a set of duplicates:
--
-- Merged or trimmed:  We compare the leftmost coordinates and the
-- aligned length.  If the library prep is strand-preserving, we also
-- compare the strand.
--
-- Paired: We compare both left-most coordinates (b_pos and b_mpos).  If
-- the library prep is strand-preserving, only first-mates can be
-- duplicates of first-mates.  Else a first-mate can be the duplicate of
-- a second-mate.  There may be pairs with one unmapped mate.  This is
-- not a problem as they get assigned synthetic coordinates and will be
-- handled smoothly.
--
-- True singles:  We compare only the leftmost coordinate.  It does not
-- matter if the library prep is strand-preserving, the strand always
-- matters.
--
-- Across these classes, we can see more duplicates:
--
-- Merged/trimmed and paired:  these can be duplicates if the merging
-- failed for the pair.  We would need to compare the outer coordinates
-- of the merged reads to the two 5' coordinates of the pair.  However,
-- since we don't have access to the mate, we cannot actually do
-- anything right here.  This case should be solved externally by
-- merging those pairs that overlap in coordinate space.
--
-- Single and paired:  in the single case, we only have one coordinate
-- to compare.  This will inevitably lead to trouble, as we could find
-- that the single might be the duplicate of two pairs, but those two
-- pairs are definitely not duplicates of each other.  We solve it by
-- removing the single read(s).
--
-- Single and merged/trimmed:  same trouble as in the single+paired
-- case.  We remove the single to solve it.
--
--
-- In principle, we might want to allow some wiggle room in the
-- coordinates.  So far, this has not been implemented.  It adds the
-- complication that groups of separated reads can turn into a set of
-- duplicates because of the appearance of a new reads.  Needs some
-- thinking about... or maybe it's not too important.
--
-- Once a set of duplicates is collected, we perform a majority vote on
-- the correct CIGAR line.  Of all those reads that agree on this CIGAR
-- line, a consensus is called, quality scores are adjusted and clamped
-- to a maximum, the MD field is updated and the XP field is assigned
-- the number of reads in the original cluster.  The new MAPQ becomes
-- the RMSQ of the map qualities of all reads.
--
-- Treatment of Read Groups:  We generalize by providing a "label"
-- function; only reads that have the same label are considered
-- duplicates of each other.  The typical label function would extract
-- read groups, libraries or samples.

rmdup :: (Monad m, Ord l) => (BamRaw -> l) -> Bool -> Collapse -> Enumeratee [BamRaw] [BamRaw] m r
rmdup label strand_preserved collapse_cfg =
    -- Easiest way to go about this:  We simply collect everything that
    -- starts at some specific coordinate and group it appropriately.
    -- Treat the groups separately, output, go on.
    check_sort ><> mapGroups rmdup_group ><> check_sort
  where
    rmdup_group = nice_sort . do_rmdup label strand_preserved collapse_cfg
    same_pos u v = br_cpos u == br_cpos v
    br_cpos br = (br_rname br, br_pos br)

    nice_sort x = sortBy (comparing br_l_seq) x

    mapGroups f o = I.tryHead >>= maybe (return o) (\a -> eneeCheckIfDone (mg1 f a []) o)
    mg1 f a acc k = I.tryHead >>= \mb -> case mb of
                        Nothing -> return . k . Chunk . f $ a:acc
                        Just b | same_pos a b -> mg1 f a (b:acc) k
                               | otherwise -> eneeCheckIfDone (mg1 f b []) . k . Chunk . f $ a:acc

check_sort :: Monad m => Enumeratee [BamRaw] [BamRaw] m a
check_sort out = I.tryHead >>= maybe (return out) (\a -> eneeCheckIfDone (step a) out)
  where
    step a k = I.tryHead >>= maybe (return . k $ Chunk [a]) (step' a k)
    step' a k b | (br_rname a, br_pos a) > (br_rname b, br_pos b) = fail msg
                | otherwise = eneeCheckIfDone (step b) . k $ Chunk [a]
    msg = "rmdup: input must be sorted"

{- Unmapped fragments should not be considered to be duplicates of
   mapped fragments.  The "unmapped" flag can serve for that:  while
   there are two classes of "unmapped reads (those that are not mapped
   and those that are mapped to an invalid position), the two sets will
   always have different coordinates.  (Unfortunately, correct duplicate
   removal now relies on correct "unmapped" and "mate unmapped" flags,
   and we didn't have those until four hours ago...)

   . Other definitions (e.g. lack of CIGAR) don't work, because that
     information won't be available for the mate.

   . This would amount to making the "unmapped" flag part of the
     coordinate, but samtools is not going to take it into account when
     sorting.

   . Instead, both flags become part of the "mate pos" grouping criterion.

 - First Mates should (probably) not be considered duplicates of Second
   Mates.  This is unconditionally true for libraries with A/B-style
   adapters (definitely 454, probably Mathias' ds protocol) and the ss
   protocol, it is not true for fork adapter protocols (standard
   Illumina).  So it's an option for now, which was apparently turned
   off for the Anuital Man.

   . Taking the option out might simplify stuff.  Do it?

 - This code ignores read groups, but it will do a majority vote on the
   RG field and call consensi for the index sequences.  If you believe
   that duplicates across read groups are impossible, you must call it
   with an appropriately filtered stream.
-}

do_rmdup :: Ord l => (BamRaw -> l) -> Bool -> Collapse -> [BamRaw] -> [BamRaw]
do_rmdup label strand_preserved Collapse{..} =
    concatMap do_rmdup1 . M.elems . accumMap label id
  where
    do_rmdup1 rds = let (ss,r1) = true_singles'
                        (ms,r2) = merged'
                        (ps,r3) = pairs'
                        ua      = unaligned'
                    in map project $ results ss ua ms ps ++ originals
                            (leftovers ss ua ms ps ++ r1 ++ r2 ++ r3)
      where
        -- Treatment of Half-Aligned Pairs (meaning one known
        -- coordinate, the validity of the alignments is immaterial)

        -- They are to be treated as follows:  Given that only one
        -- coordinate is known (5' of the aligned mate), we want to
        -- treat them like true singles.  But the unaligned mate should
        -- be kept if possible, though it should not contribute to a
        -- consensus sequence.  We assume nothing about the unaligned
        -- mate, not even that it *shouldn't* be aligned, never mind the
        -- fact that it *couldn't* be.  (*groan*)
        --
        -- Therefore, aligned reads with unaligned mate go to the same
        -- potential duplicate set as true singletons.  If a pair exists
        -- that might be a duplicate of those, the singletons and
        -- half-aligned pairs are removed.  Else a consensus is computed
        -- and remains for the aligned mate.
        --
        -- The unaligned mates end up in the same place (therefore we
        -- see them and can treat them locally).  We cannot call a
        -- consensus (these molecules may well have different length),
        -- so we select one, doesn't really matter which, and since
        -- we're treating both mates here, it doesn't even need to
        -- be reproducible without local information.
        --
        -- So to get both:  collect aligned reads with unaligned mates
        -- together with aligned true singles.  Collect the unaligned
        -- mates, which necessarily have the exact same alignment
        -- coordinates, separately.  If we don't find a matching true
        -- pair (that case is already handled smoothly), we keep the
        -- highest quality unaligned read, combine with the consensus of
        -- the half-aligned mates and singletons, and give it the
        -- lexically smallest name of the half-aligned pairs.

        (raw_pairs, raw_singles)       = partition br_isPaired rds
        (merged, true_singles)         = partition br_isMergeTrimmed raw_singles

        (pairs, raw_half_pairs)        = partition is_totally_aligned raw_pairs
        (half_unaligned, half_aligned) = partition br_isUnmapped raw_half_pairs

        is_totally_aligned b = not (br_isUnmapped b || br_isMateUnmapped b)

        mkMap f x = let m1 = M.map collapse $ accumMap f inject x
                    in (M.map fst m1, concatMap snd $ M.elems m1)

        pairs'        = mkMap (\b -> (br_mate_pos b,   strand_preserved && br_isReversed  b
                                                     , strand_preserved && br_isFirstMate b)) pairs
        merged'       = mkMap (\b -> (br_aln_length b, strand_preserved && br_isReversed  b)) merged

        true_singles' = mkMap (\b -> (                 strand_preserved && br_isReversed  b)) (true_singles++half_aligned)


        -- What to emit for unaligned mates:  we only emit a mate if we
        -- did not call a consensus, that is, if there was exactly one
        -- half-aligned read and no true singleton.  In that case, there
        -- is also exactly one mate to use and the whole input got
        -- passed through.  Else we force the result to be a singleton
        -- (by clearing appropriate flags) and all unaligned mates are
        -- passed through and marked as duplicates.

        unaligned'    = accumMap f inject half_unaligned
            where f b = strand_preserved && br_isReversed b

        results ss ua ms ps = merge_singles ss ua ( [ (rev, v) | ((_, rev),        v) <- M.toList ms ]
                                                 ++ [ (rev, v) | ((_, rev, False), v) <- M.toList ps ] )
                                               ++ [ either id id v | ((_, _, True), v) <- M.toList ps ]

        leftovers ss ua ms ps = unmerged_singles ss ua ( [ (rev, v) | ((_, rev),        v) <- M.toList ms ]
                                                      ++ [ (rev, v) | ((_, rev, False), v) <- M.toList ps ] )

    -- At every site, we need to check if were producing a consensus or
    -- if were outputting an original read.  If the latter case, its
    -- mate should be produced as ordinary else, in the former case, the
    -- bunch of unmapped mates should be produced as duplicated reads.

    merge_singles m _ [                 ] = map (either id id) $ M.elems m
    merge_singles m n ((k,Right v) : kvs) = case M.lookup k m of
            Nothing -> make_singleton v : merge_singles m (M.delete k n) kvs
            Just  w -> add_xp_of (either id id w) v : merge_singles (M.delete k m) (M.delete k n) kvs
    merge_singles m n ((k,Left  v) : kvs) = case (M.lookup k m, M.lookup k n) of
            (Nothing, Nothing) -> v : merge_singles (M.delete k m) (M.delete k n) kvs
            (Nothing, Just vs) -> v : vs ++ merge_singles (M.delete k m) (M.delete k n) kvs
            (Just  w,       _) -> add_xp_of (either id id w) v : merge_singles (M.delete k m) (M.delete k n) kvs

    unmerged_singles _ n [           ] = concat $ M.elems n
    -- the singleton is an actual consensus: emit *all* unmapped mates
    unmerged_singles m n ((k,Right _) : kvs) = case M.lookup k n of
            Nothing -> unmerged_singles (M.delete k m) (M.delete k n) kvs
            Just vs -> vs ++ unmerged_singles (M.delete k m) (M.delete k n) kvs
    -- the singleton was passed through *and* there is no pair we are
    -- going to use instead: do not pass the unmapped mates (there can
    -- be only one, and it was passed in merge_singles)
    unmerged_singles m n ((k,Left  _) : kvs) = case (M.lookup k m, M.findWithDefault [] k n) of
            (Nothing,  _) -> unmerged_singles (M.delete k m) (M.delete k n) kvs
            (Just  w, vs) -> (either id id w) : vs ++ unmerged_singles (M.delete k m) (M.delete k n) kvs

    br_mate_pos       br = (br_mrnm br, br_mpos br, br_isUnmapped br, br_isMateUnmapped br)
    br_isMergeTrimmed br = let ef = br_extAsInt (br_extAsInt 0 "XF" br) "FF" br
                           in (ef .&. flagTrimmed .|. flagMerged) /= 0

accumMap :: Ord k => (a -> k) -> (a -> v) -> [a] -> M.Map k [v]
accumMap f g = go M.empty
  where
    go m [    ] = m
    go m (a:as) = let ws = M.findWithDefault [] (f a) m
                  in ws `seq` go (M.insert (f a) (g a:ws) m) as


{- We need to deal sensibly with each field, but different fields have
   different needs.  We can take the value from the first read to
   preserve determinism or because all reads should be equal anyway,
   aggregate over all reads computing either RMSQ or the most common
   value, delete a field because it wouldn't make sense anymore or
   because doing something sensible would be hard and we're going to
   ignore it anyway, or we calculate some special value; see below.
   Unknown fields will be taken from the first read, which seems to be a
   safe default.

   QNAME and most fields              taken from first
   FLAG qc fail                       majority vote
        dup                           deleted
   MAPQ                               rmsq
   CIGAR, SEQ, QUAL, MD, NM, XP       generated
   XA                                 concatenate all
   XI/YI, XJ/YJ                       compute consensus

   BQ, CM, FZ, Q2, R2, XM, XO, XG, YQ, EN
         deleted because they would become wrong

   CQ, CS, E2, FS, OQ, OP, OC, U2, H0, H1, H2, HI, NH, IH, ZQ
         delete because they will be ignored anyway

   AM, AS, MQ, PQ, SM, UQ
         compute rmsq

   X0, X1, XT, XS, XF, XE, BC, LB, RG
         majority vote -}

addXPOf :: BamRec -> BamRec -> BamRec
addXPOf w v = v { b_exts = M.insert "XP" (Int $ extAsInt 1 "XP" w `oplus` extAsInt 1 "XP" v) (b_exts v) }

addXPOf' :: BamRaw -> BamRaw -> BamRaw
addXPOf' w v = replaceXP (br_extAsInt 1 "XP" w `oplus` br_extAsInt 1 "XP" v) v

do_collapse :: Word8 -> [BamRec] -> (Either BamRec BamRec, [BamRec])
do_collapse maxq [br] | B.all (<= maxq) (b_qual br) = ( Left     br, [  ] )     -- no modifcation, pass through
                      | otherwise                   = ( Right lq_br, [br] )     -- qualities reduces, must keep original
  where
    lq_br = br { b_qual  = B.map (min maxq) $ b_qual br
               , b_virtual_offset = 0
               , b_qname = b_qname br `B.snoc` 99 }

do_collapse maxq  brs = ( Right b0 { b_exts  = modify_extensions $ b_exts b0
                                   , b_flag  = failflag .&. complement flagDuplicate
                                   , b_mapq  = rmsq $ map b_mapq brs'
                                   , b_cigar = Cigar cigar'
                                   , b_seq   = V.fromList $ cons_seq
                                   , b_qual  = B.pack cons_qual
                                   , b_qname = b_qname b0 `B.snoc` 99
                                   , b_virtual_offset = 0 }, brs )              -- many modifications, must keep everything
  where
    b0 = minimumBy (comparing b_qname) brs
    most_fail = 2 * length (filter isFailsQC brs) > length brs
    failflag | most_fail = b_flag b0 .|. flagFailsQC
             | otherwise = b_flag b0 .&. complement flagFailsQC

    rmsq xs = round $ sqrt (x::Double)
      where
        x = fromIntegral (sum (map (\y->y*y) xs)) / genericLength xs

    maj vs = head . maximumBy (comparing length) . group . sort $ vs
    nub' = concatMap head . group . sort

    -- majority vote on the cigar lines, then filter
    cigar' = maj $ map (unCigar . b_cigar) brs
    brs' = filter ((==) cigar' . unCigar . b_cigar) brs

    (cons_seq, cons_qual) = unzip $ map (consensus maxq) $ transpose $ map to_pairs brs'

    add_index k1 k2 | null inputs = id
                    | otherwise = M.insert k1 (Text $ T.pack $ show conss) .
                                  M.insert k2 (Text $ B.pack $ map (+33) consq)
      where
        inputs = [ zip (map toNucleotide $ T.unpack sq) qs
                 | es <- map b_exts brs
                 , Text sq <- maybe [] (:[]) $ M.lookup k1 es
                 , let qs = case M.lookup k2 es of
                                Just (Text t) -> map (subtract 33) $ B.unpack t
                                _             -> repeat 23 ]
        (conss,consq) = unzip $ map (consensus 93) $ transpose $ inputs


    to_pairs b | B.null (b_qual b) = zip (V.toList $ b_seq b) (repeat 23)   -- error rate of ~0.5%
               | otherwise         = zip (V.toList $ b_seq b) (B.unpack $ b_qual b)

    md' = case [ (b_seq b,md) | b <- brs', Just md <- [ getMd b ] ] of
                [             ] -> []
                (seq1, md1) : _ -> mk_new_md cigar' md1 (V.toList seq1) cons_seq
    nm' = sum $ [ n | (Ins,n) <- cigar' ] ++ [ n | (Del,n) <- cigar' ] ++ [ 1 | MdRep _ <- md' ]
    xa' = nub' [ T.split ';' xas | Just (Text xas) <- map (M.lookup "XA" . b_exts) brs ]

    modify_extensions es = foldr ($!) es $
        [ let vs = [ v | Just v <- map (M.lookup k . b_exts) brs ]
          in if null vs then id else M.insert k $! maj vs | k <- do_maj ] ++
        [ let vs = [ v | Just (Int v) <- map (M.lookup k . b_exts) brs ]
          in if null vs then id else M.insert k $! Int (rmsq vs) | k <- do_rmsq ] ++
        [ M.delete k | k <- useless ] ++
        [ M.insert "NM" $! Int nm'
        , M.insert "XP" $! Int (foldl' (\a b -> a `oplus` extAsInt 1 "XP" b) 0 brs)
        , if null xa' then id else M.insert "XA" $! (Text $ T.intercalate (T.singleton ';') xa')
        , if null md' then id else M.insert "MD" $! (Text $ showMd md')
        , add_index "XI" "YI"
        , add_index "XJ" "YJ" ]

    useless = words "BQ CM FZ Q2 R2 XM XO XG YQ EN CQ CS E2 FS OQ OP OC U2 H0 H1 H2 HI NH IH ZQ"
    do_rmsq = words "AM AS MQ PQ SM UQ"
    do_maj  = words "X0 X1 XT XS XF XE BC LB RG"

minViewBy :: (a -> a -> Ordering) -> [a] -> (a,[a])
minViewBy  _  [    ] = error "minViewBy on empty list"
minViewBy cmp (x:xs) = go x [] xs
  where
    go m acc [    ] = (m,acc)
    go m acc (a:as) = case m `cmp` a of GT -> go a (m:acc) as
                                        _  -> go m (a:acc) as

mk_new_md :: [(CigOp, Int)] -> [MdOp] -> [Nucleotide] -> [Nucleotide] -> [MdOp]
mk_new_md [] [] [] [] = []

mk_new_md (( _ , 0):cigs) mds osq nsq = mk_new_md cigs mds osq nsq
mk_new_md cigs (MdNum  0 : mds) osq nsq = mk_new_md cigs mds osq nsq
mk_new_md cigs (MdDel [] : mds) osq nsq = mk_new_md cigs mds osq nsq

mk_new_md ((Mat, u):cigs) (MdRep b : mds) (_:osq) (n:nsq)
    | b == n    = MdNum 1 : mk_new_md ((Mat, u-1):cigs) mds osq nsq
    | otherwise = MdRep b : mk_new_md ((Mat, u-1):cigs) mds osq nsq

mk_new_md ((Mat, u):cigs) (MdNum v : mds) (o:osq) (n:nsq)
    | o == n    = MdNum 1 : mk_new_md ((Mat, u-1):cigs) (MdNum (v-1) : mds) osq nsq
    | otherwise = MdRep o : mk_new_md ((Mat, u-1):cigs) (MdNum (v-1) : mds) osq nsq

mk_new_md ((Del, n):cigs) (MdDel bs : mds) osq nsq | n == length bs = MdDel bs : mk_new_md cigs mds osq nsq

mk_new_md ((Ins, n):cigs) md osq nsq = mk_new_md cigs md (drop n osq) (drop n nsq)
mk_new_md ((SMa, n):cigs) md osq nsq = mk_new_md cigs md (drop n osq) (drop n nsq)
mk_new_md ((HMa, _):cigs) md osq nsq = mk_new_md cigs md         osq          nsq
mk_new_md ((Pad, _):cigs) md osq nsq = mk_new_md cigs md         osq          nsq
mk_new_md ((Nop, _):cigs) md osq nsq = mk_new_md cigs md         osq          nsq

mk_new_md cigs ms osq nsq = error $ unlines
    [ "Broken MD field when trying to construct new MD!"
    , "CIGAR: " ++ show cigs
    , "MD: " ++ show ms
    , "refseq: " ++ show osq
    , "readseq: " ++ show nsq ]


consensus :: Word8 -> [ (Nucleotide, Word8) ] -> (Nucleotide, Word8)
consensus maxq nqs = if qr > 3 then (n0, qr) else (nucN,0)
  where
    accs :: UArray Nucleotide Int
    accs = accumArray (+) 0 (minBound,maxBound) [ (n,fromIntegral q) | (n,q) <- nqs ]

    (n0,q0) : (_,q1) : _ = sortBy (flip $ comparing snd) $ assocs accs
    qr = fromIntegral $ (q0-q1) `min` fromIntegral maxq


-- Cheap version: simply takes the lexically first record, adds XP field
do_cheap_collapse :: [BamRaw] -> ( Either BamRaw BamRaw, [BamRaw] )
do_cheap_collapse [b] = ( Left                     b, [] )
do_cheap_collapse  bs = ( Left $ replaceXP new_xp b0, bx )
  where
    (b0, bx) = minViewBy (comparing br_qname) bs
    new_xp   = foldl' (\a b -> a `oplus` br_extAsInt 1 "XP" b) 0 bs

replaceXP :: Int -> BamRaw -> BamRaw
replaceXP new_xp b0 = bamRaw 0 . xpcode . raw_data . mutateBamRaw b0 $ removeExt "XP"
  where
    xpcode r | new_xp == 1 = r
             | -0x80 <= new_xp && new_xp < 0 = r `B.append` B.pack [ c2w 'X', c2w 'P', c2w 'c',
                                                                     fromIntegral $ new_xp .&. 0xff ]
             | new_xp < 0x100                = r `B.append` B.pack [ c2w 'X', c2w 'P', c2w 'C',
                                                                     fromIntegral $ new_xp .&. 0xff ]
             | new_xp < 0x10000              = r `B.append` B.pack [ c2w 'X', c2w 'P', c2w 'S',
                                                                     fromIntegral $ (new_xp `shiftR`  0).&. 0xff,
                                                                     fromIntegral $ (new_xp `shiftR`  8) .&. 0xff ]
             | otherwise   = r `B.append` B.pack [ c2w 'X', c2w 'P', c2w 'i',
                                                   fromIntegral $ (new_xp `shiftR`  0) .&. 0xff,
                                                   fromIntegral $ (new_xp `shiftR`  8) .&. 0xff,
                                                   fromIntegral $ (new_xp `shiftR` 16) .&. 0xff,
                                                   fromIntegral $ (new_xp `shiftR` 24) .&. 0xff ]

oplus :: Int -> Int -> Int
_ `oplus` (-1) = -1
(-1) `oplus` _ = -1
a `oplus` b = a + b

