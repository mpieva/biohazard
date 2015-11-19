{-# LANGUAGE ExistentialQuantification, RecordWildCards, NamedFieldPuns, BangPatterns, FlexibleContexts #-}
module Bio.Bam.Rmdup(
            rmdup, Collapse, cons_collapse, cheap_collapse,
            cons_collapse_keep, cheap_collapse_keep,
            check_sort, normalizeTo, wrapTo
    ) where

import Bio.Bam.Evan                     ( removeWarts )
import Bio.Bam.Header
import Bio.Bam.Raw
import Bio.Bam.Rec
import Bio.Base
import Bio.Iteratee
import Control.Monad                    ( when )
import Data.Bits
import Data.List
import Data.Ord                         ( comparing )

import qualified Data.ByteString        as B
import qualified Data.ByteString.Unsafe as B
import qualified Data.ByteString.Char8  as T
import qualified Data.Iteratee          as I
import qualified Data.Map               as M
import qualified Data.Vector.Generic    as V
import qualified Data.Vector.Unboxed    as U

-- | Uniform treatment of raw and parsed BAM records.  Might grow into
-- something even bigger.

class BAMREC a where
    inject          :: BamRaw -> a              -- convert from BamRaw
    project         :: a -> BamRaw              -- convert to BamRaw
    is_mate_of      :: a -> a -> Bool           -- check if two records form a mate
    make_singleton  :: a -> a                   -- remove all PE related flags
    flag_dup        :: a -> a                   -- flag as duplicate
    add_xp_of       :: a -> a -> a              -- add XP field of forst read to that of second

data Collapse = forall a . BAMREC a => Collapse {
                    collapse :: [a] -> (Politics a,[a]),    -- cluster to consensus and stuff or representative and stuff
                    originals :: [a] -> [a] }               -- treatment of the redundant original reads

data Politics a = Consensus a | Representative a

fromPolitics :: Politics a -> a
fromPolitics (Consensus      a) = a
fromPolitics (Representative a) = a


cons_collapse :: Qual -> Collapse
cons_collapse maxq = Collapse (do_collapse maxq) (const [])

cons_collapse_keep :: Qual -> Collapse
cons_collapse_keep maxq = Collapse (do_collapse maxq) (map flag_dup)

cheap_collapse :: Collapse
cheap_collapse = Collapse do_cheap_collapse (const [])

cheap_collapse_keep :: Collapse
cheap_collapse_keep = Collapse do_cheap_collapse (map flag_dup)

instance BAMREC BamRaw where
    inject  = id
    project = id
    make_singleton b = mutateBamRaw b $ setFlag (br_flag b .&. complement pflags)
    flag_dup       b = mutateBamRaw b $ setFlag (br_flag b .|. flagDuplicate)
    add_xp_of    w v = replaceXP (br_extAsInt 1 "XP" w `oplus` br_extAsInt 1 "XP" v) v
    is_mate_of   a b = br_qname a == br_qname b && br_isPaired a && br_isPaired b && br_isFirstMate a == br_isSecondMate b

instance BAMREC BamRec where
    inject  = removeWarts . decodeBamEntry
    project = encodeBamEntry
    make_singleton b = b { b_flag = b_flag b .&. complement pflags }
    flag_dup       b = b { b_flag = b_flag b .|. flagDuplicate }
    add_xp_of    w v = v { b_exts = updateE "XP" (Int $ extAsInt 1 "XP" w `oplus` extAsInt 1 "XP" v) (b_exts v) }
    is_mate_of   a b = b_qname a == b_qname b && isPaired a && isPaired b && isFirstMate a == isSecondMate b

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
    check_sort "input must be sorted for rmdup to work" ><>
    mapGroups rmdup_group ><>
    check_sort "internal error, output isn't sorted anymore"
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

check_sort :: Monad m => String -> Enumeratee [BamRaw] [BamRaw] m a
check_sort msg out = I.tryHead >>= maybe (return out) (\a -> eneeCheckIfDone (step a) out)
  where
    step a k = I.tryHead >>= maybe (return . k $ Chunk [a]) (step' a k)
    step' a k b | (br_rname a, br_pos a) > (br_rname b, br_pos b) = fail $ "rmdup: " ++ msg
                | otherwise = eneeCheckIfDone (step b) . k $ Chunk [a]


{- | Workhorse for duplicate removal.

 - Unmapped fragments should not be considered to be duplicates of
   mapped fragments.  The /unmapped/ flag can serve for that:  while
   there are two classes of /unmapped/ reads (those that are not mapped
   and those that are mapped to an invalid position), the two sets will
   always have different coordinates.  (Unfortunately, correct duplicate
   removal now relies on correct /unmapped/ and /mate unmapped/ flags,
   and we don't get them from unmodified BWA.  So correct operation
   requires patched BWA or a run of @bam-fixpair@.)

   (1) Other definitions (e.g. lack of CIGAR) don't work, because that
       information won't be available for the mate.

   (2) This would amount to making the /unmapped/ flag part of the
       coordinate, but samtools is not going to take it into account
       when sorting.

   (3) Instead, both flags become part of the /mate pos/ grouping
       criterion.

 - First Mates should (probably) not be considered duplicates of Second
   Mates.  This is unconditionally true for libraries with A\/B-style
   adapters (definitely 454, probably Mathias' ds protocol) and the ss
   protocol, it is not true for fork adapter protocols (vanilla Illumina
   protocol).  So it has to be an option, which would ideally be derived
   from header information.

 - This code ignores read groups, but it will do a majority vote on the
   @RG@ field and call consensi for the index sequences.  If you believe
   that duplicates across read groups are impossible, you must call it
   with an appropriately filtered stream.

 - Half-Aligned Pairs (meaning one known coordinate, while the validity
   of the alignments is immaterial) are rather complicated:

   (1) Given that only one coordinate is known (5' of the aligned mate),
       we want to treat them like true singles.  But the unaligned mate
       should be kept if possible, though it should not contribute to a
       consensus sequence.  We assume nothing about the unaligned mate,
       not even that it /shouldn't/ be aligned, never mind the fact that
       it /couldn't/ be.  (The difference is in the finite abilities of
       real world aligners, naturally.)

   (2) Therefore, aligned reads with unaligned mates go to the same
       potential duplicate set as true singletons.  If at least one pair
       exists that might be a duplicate of those, all singletons and
       half-aligned mates are discarded.  Else a consensus is computed
       and replaces the aligned mates.

   (3) The unaligned mates end up in the same place in a BAM stream as
       the aligned mates (therefore we see them and can treat them
       locally).  We cannot call a consensus, since these molecules may
       well have different length, so we select one.  It doesn't really
       matter which one is selected, and since we're treating both mates
       at the same time, it doesn't even need to be reproducible without
       local information.  This is made to be the mate of the consensus.

   (4) See 'merge_singles' for how it's actually done.
-}

do_rmdup :: Ord l => (BamRaw -> l) -> Bool -> Collapse -> [BamRaw] -> [BamRaw]
do_rmdup label strand_preserved Collapse{..} =
    concatMap do_rmdup1 . M.elems . accumMap label id
  where
    do_rmdup1 rds = map project $ results ++ originals (leftovers ++ r1 ++ r2 ++ r3)
      where
        (results, leftovers) = merge_singles singles' unaligned' $
                [ (str, fromPolitics br) | ((_,str  ),br) <- M.toList merged' ] ++
                [ (str, fromPolitics br) | ((_,str,_),br) <- M.toList pairs' ]

        (raw_pairs, raw_singles)       = partition br_isPaired rds
        (merged, true_singles)         = partition br_isMergeTrimmed raw_singles

        (pairs, raw_half_pairs)        = partition br_totally_aligned raw_pairs
        (half_unaligned, half_aligned) = partition br_isUnmapped raw_half_pairs

        mkMap f x = let m1 = M.map collapse $ accumMap f inject x
                    in (M.map fst m1, concatMap snd $ M.elems m1)

        (pairs',r1)   = mkMap (\b -> (br_mate_pos b,   br_strand b, br_mate b)) pairs
        (merged',r2)  = mkMap (\b -> (br_aln_length b, br_strand b))            merged
        (singles',r3) = mkMap                          br_strand (true_singles++half_aligned)
        unaligned'    = accumMap br_strand inject half_unaligned

        br_strand b = strand_preserved && br_isReversed   b
        br_mate   b = strand_preserved && br_isFirstMate  b


-- | Merging information about true singles, merged singles,
-- half-aligned pairs, actually aligned pairs.
--
-- We collected aligned reads with unaligned mates together with aligned
-- true singles (@singles@).  We collected the unaligned mates, which
-- necessarily have the exact same alignment coordinates, separately
-- (@unaligned@).  If we don't find a matching true pair (that case is
-- already handled smoothly), we keep the highest quality unaligned
-- mate, pair it with the consensus of the aligned mates and aligned
-- singletons, and give it the lexically smallest name of the
-- half-aligned pairs.

-- NOTE:  I need to decide when to run 'make_singleton'.  Basically,
-- when we call a consensus for half-aligned pairs and keep
-- everything(?).  Then we don't have a mate for the consensus... though
-- we could decide to duplicate one mate read to get it.

merge_singles :: BAMREC a
              => M.Map Bool (Politics a)                -- strand --> true singles & half aligned
              -> M.Map Bool [a]                         -- strand --> half unaligned
              -> [ (Bool, a) ]                          -- strand --> paireds & mergeds
              -> ([a],[a])                              -- results, leftovers

merge_singles singles unaligneds = go
  where
    -- Say we generated a consensus or passed something through.  If
    -- there is a singleton consensus with the same strand, we should
    -- add in its XP field and discard it.  If there is a singleton
    -- representative, we add in its XP field and put it into the
    -- leftovers.  If there is unaligned stuff here that has the same
    -- strand, it goes to the leftovers.
    go ( (str, v) : paireds) =
        let (r,l) = merge_singles (M.delete str singles) (M.delete str unaligneds) paireds
            unal  = M.findWithDefault [] str unaligneds ++ l

        in case M.lookup str singles of
            Nothing                 -> (             v : r,     unal )
            Just (Consensus      w) -> ( add_xp_of w v : r,     unal )      -- XXX do we need this w?!
            Just (Representative w) -> ( add_xp_of w v : r, w : unal )

    -- No more pairs, delegate the problem
    go [] = merge_halves unaligneds (M.toList singles)


-- | Merging of half-aligned reads.  The first argument is a map of
-- unaligned reads (their mates are aligned to the current position),
-- the second is a list of reads that are aligned (their mates are not
-- aligned).
--
-- So, suppose we're looking at a 'Representative' that was passed
-- through.  We need to emit it along with its mate, which may be hidden
-- inside a list.  (Alternatively, we could force it to single, but that
-- fails if we're passing everything along somehow.)
--
-- Suppose we're looking at a 'Consensus'.  We could pair it with some
-- mate (which we'd need to duplicate), or we could turn it into a
-- singleton.  Duplication is ugly, so in this case, we force it to
-- singleton.

merge_halves :: BAMREC a
             => M.Map Bool [a]                          -- strand --> half unaligned
             -> [(Bool, Politics a)]                    -- strand --> true singles & half aligned
             -> ([a],[a])                               -- results, leftovers

-- Emitting a consensus: make it a single.  Nothing goes to leftovers;
-- we may still need it for something else to be emitted.  (While that
-- would be strange, making sure the BAM file stays completely valid is
-- probably better.)
merge_halves unaligneds ((_, Consensus v) : singles) =
    let (r,l) = merge_halves unaligneds singles
    in (make_singleton v : r, l)

-- Emitting a representative:  find the mate in the list of unaligned
-- reads (take up to one match to be robust), and emit that, too, as a
-- result.  Everything else goes to leftovers.  If the representative
-- happens to be unpaired, no mate is found and that case therefore is
-- handled smoothly.
merge_halves unaligneds ((str, Representative v) : singles) =
    let (r,l) = merge_halves (M.delete str unaligneds) singles
        (same,diff) = partition (is_mate_of v) $ M.findWithDefault [] str unaligneds
    in (v : take 1 same ++ r, drop 1 same ++ diff ++ l)

-- No more singles, all unaligneds are leftovers.
merge_halves unaligneds [] = ( [], concat $ M.elems unaligneds )




type MPos = (Refseq, Int, Bool, Bool)

br_mate_pos :: BamRaw -> MPos
br_mate_pos br = (br_mrnm br, br_mpos br, br_isUnmapped br, br_isMateUnmapped br)

br_totally_aligned :: BamRaw -> Bool
br_totally_aligned br = not (br_isUnmapped br || br_isMateUnmapped br)


accumMap :: Ord k => (a -> k) -> (a -> v) -> [a] -> M.Map k [v]
accumMap f g = go M.empty
  where
    go m [    ] = m
    go m (a:as) = let ws = M.findWithDefault [] (f a) m ; g' = g a
                  in g' `seq` go (M.insert (f a) (g':ws) m) as


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

   BQ, CM, FZ, Q2, R2, XM, XO, XG, YQ, EN
         deleted because they would become wrong

   CQ, CS, E2, FS, OQ, OP, OC, U2, H0, H1, H2, HI, NH, IH, ZQ
         delete because they will be ignored anyway

   AM, AS, MQ, PQ, SM, UQ
         compute rmsq

   X0, X1, XT, XS, XF, XE, BC, LB, RG, XI, YI, XJ, YJ
         majority vote -}

do_collapse :: Qual -> [BamRec] -> (Politics BamRec, [BamRec])
do_collapse (Q maxq) [br] | B.all (<= maxq) (b_qual br) = ( Representative br, [  ] )     -- no modifcation, pass through
                          | otherwise                   = ( Consensus   lq_br, [br] )     -- qualities reduced, must keep original
  where
    lq_br = br { b_qual  = B.map (min maxq) $ b_qual br
               , b_virtual_offset = 0
               , b_qname = b_qname br `B.snoc` c2w 'c' }

do_collapse maxq  brs = ( Consensus b0 { b_exts  = modify_extensions $ b_exts b0
                                       , b_flag  = failflag .&. complement flagDuplicate
                                       , b_mapq  = Q $ rmsq $ map (unQ . b_mapq) $ good brs
                                       , b_cigar = Cigar cigar'
                                       , b_seq   = V.fromList $ map fst cons_seq_qual
                                       , b_qual  = B.pack $ map (unQ . snd) cons_seq_qual
                                       , b_qname = b_qname b0 `B.snoc` 99
                                       , b_virtual_offset = 0 }, brs )              -- many modifications, must keep everything
  where
    !b0 = minimumBy (comparing b_qname) brs
    !most_fail = 2 * length (filter isFailsQC brs) > length brs
    !failflag | most_fail = b_flag b0 .|. flagFailsQC
              | otherwise = b_flag b0 .&. complement flagFailsQC

    rmsq xs = case foldl' (\(!n,!d) x -> (n + fromIntegral x * fromIntegral x, d + 1)) (0,0) xs of
        (!n,!d) -> round $ sqrt $ (n::Double) / fromIntegral (d::Int)

    maj xs = head . maximumBy (comparing length) . group . sort $ xs
    nub' = concatMap head . group . sort

    -- majority vote on the cigar lines, then filter
    !cigar' = maj $ map (unCigar . b_cigar) brs
    good = filter ((==) cigar' . unCigar . b_cigar)

    cons_seq_qual = [ consensus maxq [ (V.unsafeIndex (b_seq b) i, Q q)
                                     | b <- good brs, let q = if B.null (b_qual b) then 23 else B.unsafeIndex (b_qual b) i ]
                    | i <- [0 .. len - 1] ]
        where !len = V.length . b_seq . head $ good brs

    md' = case [ (b_seq b,md,b) | b <- good brs, Just md <- [ getMd b ] ] of
                [               ] -> []
                (seq1, md1,b) : _ -> case mk_new_md cigar' md1 (V.toList seq1) (map fst cons_seq_qual) of
                    Right x -> x
                    Left (MdFail cigs ms osq nsq) -> error $ unlines
                                    [ "Broken MD field when trying to construct new MD!"
                                    , "QNAME: " ++ show (b_qname b)
                                    , "POS:   " ++ shows (unRefseq (b_rname b)) ":" ++ show (b_pos b)
                                    , "CIGAR: " ++ show cigs
                                    , "MD:    " ++ show ms
                                    , "refseq:  " ++ show osq
                                    , "readseq: " ++ show nsq ]


    nm' = sum $ [ n | (Ins,n) <- cigar' ] ++ [ n | (Del,n) <- cigar' ] ++ [ 1 | MdRep _ <- md' ]
    xa' = nub' [ T.split ';' xas | Just (Text xas) <- map (lookup "XA" . b_exts) brs ]

    modify_extensions es = foldr ($!) es $
        [ let vs = [ v | Just v <- map (lookup k . b_exts) brs ]
          in if null vs then id else updateE k $! maj vs | k <- do_maj ] ++
        [ let vs = [ v | Just (Int v) <- map (lookup k . b_exts) brs ]
          in if null vs then id else updateE k $! Int (rmsq vs) | k <- do_rmsq ] ++
        [ deleteE k | k <- useless ] ++
        [ updateE "NM" $! Int nm'
        , updateE "XP" $! Int (foldl' (\a b -> a `oplus` extAsInt 1 "XP" b) 0 brs)
        , if null xa' then id else updateE "XA" $! (Text $ T.intercalate (T.singleton ';') xa')
        , if null md' then id else updateE "MD" $! (Text $ showMd md') ]

    useless = words "BQ CM FZ Q2 R2 XM XO XG YQ EN CQ CS E2 FS OQ OP OC U2 H0 H1 H2 HI NH IH ZQ"
    do_rmsq = words "AM AS MQ PQ SM UQ"
    do_maj  = words "X0 X1 XT XS XF XE BC LB RG XI XJ YI YJ"

minViewBy :: (a -> a -> Ordering) -> [a] -> (a,[a])
minViewBy  _  [    ] = error "minViewBy on empty list"
minViewBy cmp (x:xs) = go x [] xs
  where
    go m acc [    ] = (m,acc)
    go m acc (a:as) = case m `cmp` a of GT -> go a (m:acc) as
                                        _  -> go m (a:acc) as

data MdFail = MdFail [(CigOp, Int)] [MdOp] [Nucleotides] [Nucleotides]

mk_new_md :: [(CigOp, Int)] -> [MdOp] -> [Nucleotides] -> [Nucleotides] -> Either MdFail [MdOp]
mk_new_md = mk_new_md' []

mk_new_md' :: [MdOp] -> [(CigOp, Int)] -> [MdOp] -> [Nucleotides] -> [Nucleotides] -> Either MdFail [MdOp]
mk_new_md' acc [] [] [] [] = Right $ normalize [] acc
    where
        normalize          a  (MdNum  0:os) = normalize               a  os
        normalize (MdNum n:a) (MdNum  m:os) = normalize (MdNum  (n+m):a) os
        normalize          a  (MdDel []:os) = normalize               a  os
        normalize (MdDel u:a) (MdDel  v:os) = normalize (MdDel (v++u):a) os
        normalize          a  (       o:os) = normalize            (o:a) os
        normalize          a  [           ] = a

mk_new_md' acc (( _ , 0):cigs)  mds  osq nsq = mk_new_md' acc cigs mds osq nsq
mk_new_md' acc cigs (MdNum  0 : mds) osq nsq = mk_new_md' acc cigs mds osq nsq
mk_new_md' acc cigs (MdDel [] : mds) osq nsq = mk_new_md' acc cigs mds osq nsq

mk_new_md' acc ((Mat, u):cigs) (MdRep b : mds) (_:osq) (n:nsq)
    | b == n    = mk_new_md' (MdNum 1 : acc) ((Mat, u-1):cigs) mds osq nsq
    | otherwise = mk_new_md' (MdRep b : acc) ((Mat, u-1):cigs) mds osq nsq

mk_new_md' acc ((Mat, u):cigs) (MdNum v : mds) (o:osq) (n:nsq)
    | o == n    = mk_new_md' (MdNum 1 : acc) ((Mat, u-1):cigs) (MdNum (v-1) : mds) osq nsq
    | otherwise = mk_new_md' (MdRep o : acc) ((Mat, u-1):cigs) (MdNum (v-1) : mds) osq nsq

mk_new_md' acc ((Del, n):cigs) (MdDel bs : mds) osq nsq | n == length bs = mk_new_md' (MdDel bs : acc)    cigs               mds  osq nsq
mk_new_md' acc ((Del, n):cigs) (MdDel (b:bs) : mds) osq nsq = mk_new_md' (MdDel     [b] : acc) ((Del,n-1):cigs) (MdDel    bs:mds) osq nsq
mk_new_md' acc ((Del, n):cigs) (MdRep   b    : mds) osq nsq = mk_new_md' (MdDel     [b] : acc) ((Del,n-1):cigs)              mds  osq nsq
mk_new_md' acc ((Del, n):cigs) (MdNum   m    : mds) osq nsq = mk_new_md' (MdDel [nucsN] : acc) ((Del,n-1):cigs) (MdNum (m-1):mds) osq nsq

mk_new_md' acc ((Ins, n):cigs) md osq nsq = mk_new_md' acc cigs md (drop n osq) (drop n nsq)
mk_new_md' acc ((SMa, n):cigs) md osq nsq = mk_new_md' acc cigs md (drop n osq) (drop n nsq)
mk_new_md' acc ((HMa, _):cigs) md osq nsq = mk_new_md' acc cigs md         osq          nsq
mk_new_md' acc ((Pad, _):cigs) md osq nsq = mk_new_md' acc cigs md         osq          nsq
mk_new_md' acc ((Nop, _):cigs) md osq nsq = mk_new_md' acc cigs md         osq          nsq

mk_new_md' _acc cigs ms osq nsq = Left $ MdFail cigs ms osq nsq

consensus :: Qual -> [ (Nucleotides, Qual) ] -> (Nucleotides, Qual)
consensus (Q maxq) nqs = if qr > 3 then (n0, Q qr) else (nucsN, Q 0)
  where
    accs :: U.Vector Int
    accs = U.accum (+) (U.replicate 16 0) [ (fromIntegral n, fromIntegral q) | (Ns n,Q q) <- nqs ]

    (n0,q0) : (_,q1) : _ = sortBy (flip $ comparing snd) $ zip [Ns 0 ..] $ U.toList accs
    qr = fromIntegral $ (q0-q1) `min` fromIntegral maxq


-- Cheap version: simply takes the lexically first record, adds XP field
do_cheap_collapse :: [BamRaw] -> ( Politics BamRaw, [BamRaw] )
do_cheap_collapse [b] = ( Representative                     b, [] )
do_cheap_collapse  bs = ( Representative $ replaceXP new_xp b0, bx )
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

-- | Normalize a read's alignment to fall into the canonical region
-- of [0..l].  Takes the name of the reference sequence and its length.
normalizeTo :: Seqid -> Int -> BamRaw -> BamRaw
normalizeTo nm l br = mutateBamRaw br $ do setPos (br_pos br `mod` l)
                                           setMpos (br_mpos br `mod` l)
                                           setBin (br_pos br `mod` l) (br_aln_length br)
                                           when dups_are_fine $ setMapq 37 >> removeExt "XA"
  where
    dups_are_fine  = all_match_XA (br_extAsString "XA" br)
    all_match_XA s = case T.split ';' s of [xa1, xa2] | T.null xa2 -> one_match_XA xa1
                                           [xa1]                   -> one_match_XA xa1
                                           _                       -> False
    one_match_XA s = case T.split ',' s of (sq:pos:_) | sq == nm   -> pos_match_XA pos ; _ -> False
    pos_match_XA s = case T.readInt s   of Just (p,z) | T.null z   -> int_match_XA p ;   _ -> False
    int_match_XA p | p >= 0    =  (p-1) `mod` l == br_pos br `mod` l && not (br_isReversed br)
                   | otherwise = (-p-1) `mod` l == br_pos br `mod` l && br_isReversed br


-- | Wraps a read to be fully contained in the canonical interval
-- [0..l].  If the read overhangs, it is duplicated and both copies are
-- suitably masked.
wrapTo :: Int -> BamRaw -> [BamRaw]
wrapTo l br = if overhangs then do_wrap else [br]
  where
    overhangs = not (br_isUnmapped br) && br_pos br < l && l < br_pos br + br_aln_length br

    do_wrap = let b = decodeBamEntry br in
              case split_ecig (l - b_pos b) $ toECig (b_cigar b) (maybe [] id $ getMd b) of
                  (left,right) -> [ encodeBamEntry $ b { b_cigar = toCigar  left }            `setMD` left
                                  , encodeBamEntry $ b { b_cigar = toCigar right, b_pos = 0 } `setMD` right ]

-- | Split an 'ECig' into two at some position.  The position is counted
-- in terms of the reference (therefore, deletions count, insertions
-- don't).  The parts that would be skipped if we were splitting lists
-- are replaced by soft masks.
split_ecig :: Int -> ECig -> (ECig, ECig)
split_ecig _    WithMD = (WithMD,       WithMD)
split_ecig _ WithoutMD = (WithoutMD, WithoutMD)
split_ecig 0       ecs = (mask_all ecs,    ecs)

split_ecig i (Ins' n ecs) = case split_ecig i ecs of (u,v) -> (Ins' n u, SMa' n v)
split_ecig i (SMa' n ecs) = case split_ecig i ecs of (u,v) -> (SMa' n u, SMa' n v)
split_ecig i (HMa' n ecs) = case split_ecig i ecs of (u,v) -> (HMa' n u, HMa' n v)
split_ecig i (Pad' n ecs) = case split_ecig i ecs of (u,v) -> (Pad' n u,        v)

split_ecig i (Mat' n ecs)
    | i >= n    = case split_ecig (i-n) ecs of (u,v) -> (Mat' n u, SMa' n v)
    | otherwise = (Mat' i $ SMa' (n-i) $ mask_all ecs, SMa' i $ Mat' (n-i) ecs)

split_ecig i (Rep' x ecs) = case split_ecig (i-1) ecs of (u,v) -> (Rep' x u, SMa' 1 v)
split_ecig i (Del' x ecs) = case split_ecig (i-1) ecs of (u,v) -> (Del' x u,        v)

split_ecig i (Nop' n ecs)
    | i >= n    = case split_ecig (i-n) ecs of (u,v) -> (Nop' n u,        v)
    | otherwise = (Nop' i $ mask_all ecs, Nop' (n-i) ecs)

mask_all :: ECig -> ECig
mask_all      WithMD = WithMD
mask_all   WithoutMD = WithoutMD
mask_all (Nop' _ ec) =          mask_all ec
mask_all (HMa' _ ec) =          mask_all ec
mask_all (Pad' _ ec) =          mask_all ec
mask_all (Del' _ ec) =          mask_all ec
mask_all (Rep' _ ec) = SMa' 1 $ mask_all ec
mask_all (Mat' n ec) = SMa' n $ mask_all ec
mask_all (Ins' n ec) = SMa' n $ mask_all ec
mask_all (SMa' n ec) = SMa' n $ mask_all ec

-- | Argh, this business with the CIGAR operations is a mess, it gets
-- worse when combined with MD.  Okay, we will support CIGAR (no "=" and
-- "X" operations) and MD.  If we have MD on input, we generate it on
-- output, too.  And in between, we break everything into /very small/
-- operations.  (Yes, the two terminating constructors are a weird
-- hack.)

data ECig = WithMD                      -- terminate, do generate MD field
          | WithoutMD                   -- terminate, don't bother with MD
          | Mat' Int ECig
          | Rep' Nucleotides ECig
          | Ins' Int ECig
          | Del' Nucleotides ECig
          | Nop' Int ECig
          | SMa' Int ECig
          | HMa' Int ECig
          | Pad' Int ECig


toECig :: Cigar -> [MdOp] -> ECig
toECig (Cigar cig) md = go cig md
  where
    go        cs  (MdNum  0:mds) = go cs mds
    go        cs  (MdDel []:mds) = go cs mds
    go ((_,0):cs)           mds  = go cs mds
    go [        ] [            ] = WithMD               -- all was fine to the very end
    go [        ]              _ = WithoutMD            -- here it wasn't fine

    go ((Mat,n):cs) (MdRep x:mds)      = Rep'   x   $ go  ((Mat,n-1):cs)             mds
    go ((Mat,n):cs) (MdNum m:mds)
       | n < m                         = Mat'   n   $ go             cs (MdNum (m-n):mds)
       | n > m                         = Mat'   m   $ go  ((Mat,n-m):cs)             mds
       | otherwise                     = Mat'   n   $ go             cs              mds
    go ((Mat,n):cs)            _       = Mat'   n   $ go'            cs

    go ((Ins,n):cs)               mds  = Ins'   n   $ go             cs              mds
    go ((Del,n):cs) (MdDel (x:xs):mds) = Del'   x   $ go  ((Del,n-1):cs)   (MdDel xs:mds)
    go ((Del,n):cs)                 _  = Del' nucsN $ go' ((Del,n-1):cs)

    go ((Nop,n):cs) mds = Nop' n $ go cs mds
    go ((SMa,n):cs) mds = SMa' n $ go cs mds
    go ((HMa,n):cs) mds = HMa' n $ go cs mds
    go ((Pad,n):cs) mds = Pad' n $ go cs mds

    -- We jump here once the MD fiels ran out early or was messed up.
    -- We no longer bother with it (this also happens if the MD isn't
    -- present to begin with).
    go' ((_,0):cs)   = go' cs
    go' [        ]   = WithoutMD                        -- we didn't have MD or it was broken

    go' ((Mat,n):cs) = Mat'   n   $ go'            cs
    go' ((Ins,n):cs) = Ins'   n   $ go'            cs
    go' ((Del,n):cs) = Del' nucsN $ go' ((Del,n-1):cs)

    go' ((Nop,n):cs) = Nop'   n   $ go' cs
    go' ((SMa,n):cs) = SMa'   n   $ go' cs
    go' ((HMa,n):cs) = HMa'   n   $ go' cs
    go' ((Pad,n):cs) = Pad'   n   $ go' cs


-- We normalize matches, deletions and soft masks, because these are the
-- operations we generate.  Everything else is either already normalized
-- or nobody really cares anyway.
toCigar :: ECig -> Cigar
toCigar = Cigar . go
  where
    go       WithMD = []
    go    WithoutMD = []

    go (Ins' n ecs) = (Ins,n) : go ecs
    go (Nop' n ecs) = (Nop,n) : go ecs
    go (HMa' n ecs) = (HMa,n) : go ecs
    go (Pad' n ecs) = (Pad,n) : go ecs
    go (SMa' n ecs) = go_sma n ecs
    go (Mat' n ecs) = go_mat n ecs
    go (Rep' _ ecs) = go_mat 1 ecs
    go (Del' _ ecs) = go_del 1 ecs

    go_sma !n (SMa' m ecs) = go_sma (n+m) ecs
    go_sma !n         ecs  = (SMa,n) : go ecs

    go_mat !n (Mat' m ecs) = go_mat (n+m) ecs
    go_mat !n (Rep' _ ecs) = go_mat (n+1) ecs
    go_mat !n         ecs  = (Mat,n) : go ecs

    go_del !n (Del' _ ecs) = go_del (n+1) ecs
    go_del !n         ecs  = (Del,n) : go ecs



-- | Create an MD field from an extended CIGAR and place it in a record.
-- We build it piecemeal (in 'go'), call out to 'addNum', 'addRep',
-- 'addDel' to make sure the operations are not generated in a
-- degenerate manner, and finally check if we're even supposed to create
-- an MD field.
setMD :: BamRec -> ECig -> BamRec
setMD b ec = case go ec of
    Just md -> b { b_exts = updateE "MD" (Text $ showMd md) (b_exts b) }
    Nothing -> b { b_exts = deleteE "MD"                    (b_exts b) }
  where
    go  WithMD      = Just []
    go  WithoutMD   = Nothing

    go (Ins' _ ecs) = go ecs
    go (Nop' _ ecs) = go ecs
    go (SMa' _ ecs) = go ecs
    go (HMa' _ ecs) = go ecs
    go (Pad' _ ecs) = go ecs
    go (Mat' n ecs) = (if n ==  0 then id else fmap (addNum n)) $ go ecs
    go (Rep' x ecs) = (if isGap x then id else fmap (addRep x)) $ go ecs
    go (Del' x ecs) = (if isGap x then id else fmap (addDel x)) $ go ecs

    addNum n (MdNum m : mds) = MdNum (n+m) : mds
    addNum n            mds  = MdNum   n   : mds

    addRep x            mds  = MdRep   x   : mds

    addDel x (MdDel y : mds) = MdDel (x:y) : mds
    addDel x            mds  = MdDel  [x]  : mds
