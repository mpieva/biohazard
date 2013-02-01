{-# LANGUAGE ExistentialQuantification #-}
module Bio.File.Bam.Rmdup(
            rmdup, Collapse,
            cons_collapse, cheap_collapse, very_cheap_collapse
    ) where

import Bio.File.Bam
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

data Collapse = forall a . Collapse (BamRaw -> Either String a) ([a] -> a) (a -> a -> a) (a -> BamRaw)

cons_collapse :: Word8 -> Collapse
cons_collapse maxq = Collapse check_flags (do_collapse maxq) addXPOf encodeBamEntry

cheap_collapse :: Collapse
cheap_collapse = Collapse check_flags do_cheap_collapse addXPOf encodeBamEntry

very_cheap_collapse :: Collapse
very_cheap_collapse = Collapse Right do_very_cheap_collapse (const id) id


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
    check_sort ><> mapGroups (either fail return . do_rmdup label strand_preserved collapse_cfg) ><> check_sort
  where
    same_pos u v = br_cpos u == br_cpos v
    br_cpos br = (br_rname br, br_pos br)

    mapGroups f o = I.tryHead >>= maybe (return o) (\a -> eneeCheckIfDone (mg1 f a []) o)
    mg1 f a acc k = I.tryHead >>= \mb -> case mb of
                        Nothing -> f (a:acc) >>= return . k . Chunk
                        Just b | same_pos a b -> mg1 f a (b:acc) k
                               | otherwise -> f (a:acc) >>= eneeCheckIfDone (mg1 f b []) . k . Chunk
    
check_sort :: Monad m => Enumeratee [BamRaw] [BamRaw] m a
check_sort out = I.tryHead >>= maybe (return out) (\a -> eneeCheckIfDone (step a) out)
  where
    step a k = I.tryHead >>= maybe (return . k $ Chunk [a]) (step' a k)
    step' a k b | (br_rname a, br_pos a) > (br_rname b, br_pos b) = fail "sorting violated"
                | otherwise = eneeCheckIfDone (step b) . k $ Chunk [a]

-- To be perfectly honest, I do not understand what these flags mean.
-- All I know it that if and when they are set, the duplicate removal
-- will go very, very wrong.
check_flags :: BamRaw -> Either String BamRec
check_flags raw | extAsInt 1 "HI" b /= 1 = Left "cannot deal with HI /= 1"
                | extAsInt 1 "IH" b /= 1 = Left "cannot deal with IH /= 1"
                | extAsInt 1 "NH" b /= 1 = Left "cannot deal with NH /= 1"
                | otherwise              = Right b
  where b = removeWarts $ decodeBamEntry raw
    

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

do_rmdup :: Ord l => (BamRaw -> l) -> Bool -> Collapse -> [BamRaw] -> Either String [BamRaw]
do_rmdup label strand_preserved (Collapse inject collapse add_xp_of finish) = 
    fmap concat . mapM do_rmdup1 . M.elems . accumMap label
  where
    do_rmdup1 rds = do ss <- true_singles'
                       ms <- merged'
                       ps <- pairs'
                       return $ map finish $ results ss ms ps
      where
        (pairs,  singles)      = partition br_isPaired rds
        (merged, true_singles) = partition br_isMergeTrimmed singles

        mkMap f = fmap (M.map collapse) . accumMapM f inject

        pairs'        = mkMap (\b -> (br_mate_pos b,   strand_preserved && br_isReversed  b 
                                                     , strand_preserved && br_isFirstMate b)) pairs
        merged'       = mkMap (\b -> (br_aln_length b, strand_preserved && br_isReversed  b)) merged

        true_singles' = mkMap (\b -> (                 strand_preserved && br_isReversed  b)) true_singles

        results ss ms ps = merge_singles ss ( [ (rev, v) | ((_, rev),        v) <- M.toList ms ]
                                           ++ [ (rev, v) | ((_, rev, False), v) <- M.toList ps ] )
                                        ++ [ v | ((_, _, True), v) <- M.toList ps ]

    merge_singles m [           ] = M.elems m
    merge_singles m ((k,v) : kvs) = case M.lookup k m of
            Nothing -> v : merge_singles m kvs
            Just  w -> add_xp_of w v : merge_singles (M.delete k m) kvs

    br_mate_pos       br = (br_mrnm br, br_mpos br, br_isUnmapped br, br_isMateUnmapped br)
    br_isMergeTrimmed br = let ef = br_extAsInt (br_extAsInt 0 "XF" br) "EF" br 
                           in (ef .&. flagTrimmed .|. flagMerged) /= 0

accumMap :: Ord b => (a -> b) -> [a] -> M.Map b [a]
accumMap f = go M.empty    
  where    
    go m [    ] = m
    go m (v:vs) = let ws = M.findWithDefault [] (f v) m 
                  in ws `seq` go (M.insert (f v) (v:ws) m) vs

accumMapM :: (Ord k, Monad m) => (a -> k) -> (a -> m v) -> [a] -> m (M.Map k [v])
accumMapM f g = go M.empty
  where    
    go m [    ] = return m
    go m (v:vs) = do v' <- g v
                     let ws = M.findWithDefault [] (f v) m 
                     ws `seq` go (M.insert (f v) (v':ws) m) vs


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

do_collapse :: Word8 -> [BamRec] -> BamRec
do_collapse maxq [br] = br { b_qual  = B.map (min maxq) $ b_qual br, b_virtual_offset = 0 }
do_collapse maxq  brs = b0 { b_exts  = modify_extensions $ b_exts b0
                           , b_flag  = failflag .&. complement flagDuplicate
                           , b_mapq  = rmsq $ map b_mapq brs'
                           , b_cigar = Cigar cigar'
                           , b_seq   = V.fromList $ cons_seq
                           , b_qual  = B.pack cons_qual
                           , b_virtual_offset = 0 }
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


-- Cheap version: simply takes the lexically first record
do_very_cheap_collapse :: [BamRaw] -> BamRaw
do_very_cheap_collapse bs = minimumBy (comparing br_qname) bs

-- Cheap version: simply takes the lexically first record, adds XP field
do_cheap_collapse :: [BamRec] -> BamRec
do_cheap_collapse bs = b0 { b_exts = new_xp $ b_exts b0 }
  where
    b0     = minimumBy (comparing b_qname) bs
    new_xp = M.insert "XP" $! Int (foldl' (\a b -> a `oplus` extAsInt 1 "XP" b) 0 bs)

oplus :: Int -> Int -> Int
_ `oplus` (-1) = -1
(-1) `oplus` _ = -1
a `oplus` b = a + b

