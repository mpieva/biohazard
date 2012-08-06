module Bio.File.Bam.Rmdup ( rmdup ) where

import Bio.File.Bam
import Bio.Iteratee

import Data.Array.Unboxed
import Data.List
import Data.Ord
import Data.Word ( Word8 )

import qualified Data.Map as M
import qualified Data.ByteString as SB
import qualified Data.Iteratee as I

-- | Removes duplicates from an aligned, sorted BAM stream.
--
-- The incoming stream must be sorted by coordinate, and we check for
-- violations of that assumption.  We cannot assume that length was
-- taken into account when sorting (samtools doesn't do so), so
-- duplicates may be separated by reads that start at the same position
-- but have different length or different strand.
--
-- We are looking at three different kinds of reads:  paired reads, true
-- single ended reads, merged or trimmed reads.  They do not overlap, by
-- definition, and are treated separately.  These conditions define a
-- set of duplicates:
--
-- Merged or trimmed:  We compare the leftmost coordinates and the
-- aligned length.  If the library prep is strand-preserving, we also
-- compare the strand.
--
-- Paired: We compare both left-most coordinates (b_pos and b_mpos).  If
-- the library prep is strand-preserving, only first-mates can be
-- duplicates of first-mates.  Else a first-mate can be the duplicate of
-- a second-mate.  In principle, pairs can be duplicates of merged
-- reads.  We do not handle this case, it should be solved externally by
-- merging those pairs that overlap in coordinate space.  Also, there
-- may be pairs with one unmapped mate.  This is not a problem as they
-- get assigned synthetic coordinates and will be handled smoothly.
--
-- True singles:  We compare only the leftmost coordinate.  It does not
-- matter if the library prep is strand-preserving.  Technically, a true
-- single could be the duplicate of a trimmed single, if the trimming is
-- unreliable.  We do not consider this possibility... some losses are
-- simply inevitable.
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

rmdup :: Monad m => Bool -> Word8 -> Enumeratee [BamRec] [BamRec] m a 
rmdup strand_preserved maxq =
    -- Easiest way to go about this:  We simply collect everything that
    -- starts at some specific coordinate and group it appropriately.
    -- Treat the groups separately, output, go on.
    check_sort ><> mapGroups (do_rmdup strand_preserved maxq) ><> check_sort
  where
    same_pos u v = b_cpos u == b_cpos v
    b_cpos br = (b_rname br, b_pos br)

    mapGroups f o = I.tryHead >>= maybe (return o) (\a -> eneeCheckIfDone (mg1 f a []) o)
    mg1 f a acc k = I.tryHead >>= \mb -> case mb of
                        Nothing -> return . k . Chunk . f $ a : acc
                        Just b | same_pos a b -> mg1 f a (b:acc) k
                               | otherwise -> eneeCheckIfDone (mg1 f b []) . k . Chunk . f $ a : acc
    
check_sort :: Monad m => Enumeratee [BamRec] [BamRec] m a
check_sort out = I.tryHead >>= maybe (return out) (\a -> eneeCheckIfDone (step a) out)
  where
    step a k = I.tryHead >>= maybe (return . k $ Chunk [a]) (step' a k)
    step' a k b | (b_rname a, b_pos a) > (b_rname b, b_pos b) = fail "sorting violated"
                | otherwise = eneeCheckIfDone (step b) . k $ Chunk [a]

do_rmdup :: Bool -> Word8 -> [BamRec] -> [BamRec]
do_rmdup strand_preserved maxq rds = map (collapse maxq) $ filter (not . null) groups
  where
    (pairs, singles) = partition isPaired rds
    (merged, true_singles) = partition isMergeTrimmed singles

    isMergeTrimmed br = case M.lookup "XF" (b_exts br) of
        Just (Int i) -> rem i 4 /= 0 ; _ -> False

    b_aln_len = cigarToAlnLen . b_cigar
    b_mate_pos br = (b_mrnm br, b_mpos br)

    group'sort f = groupBy (\a b -> f a == f b) . sortBy (comparing f)

    groups = true_singles :
             group'sort (\b -> (b_aln_len b,  strand_preserved && isReversed  b)) merged ++ 
             group'sort (\b -> (b_mate_pos b, strand_preserved && isFirstMate b)) pairs

collapse :: Word8 -> [BamRec] -> BamRec
collapse maxq [br] = br { b_qual = SB.map (min maxq) $ b_qual br, b_virtual_offset = 0 }
collapse maxq  brs = b0 { b_exts = M.insert "XP" (Int (foldl' oplus 0 (map get_xp brs))) $
                                   M.insert "MD" (Text $ showMd md') $ (b_exts b0)
                        , b_mapq = rmsq $ map b_mapq brs'
                        , b_cigar = Cigar cigar'
                        , b_seq = cons_seq
                        , b_qual = SB.pack cons_qual
                        , b_virtual_offset = 0 }
  where
    _ `oplus` (-1) = -1
    (-1) `oplus` _ = -1
    a `oplus` b = a + b

    b0 = minimumBy (comparing b_qname) brs

    rmsq xs = round $ sqrt x
      where
        x :: Double 
        x = fromIntegral (sum (map (\y->y*y) xs)) / genericLength xs

    -- majority vote on the cigar lines, then filter
    cigar' = head . maximumBy (comparing length) . group . sort $ map (unCigar . b_cigar) brs
    brs' = filter ((==) cigar' . unCigar . b_cigar) brs
    (seq1, md1) : _ = [ (b_seq b,md) | b <- brs', Just md <- [ getMd b ] ]

    (cons_seq, cons_qual) = unzip $ map (consensus maxq) $ transpose $ map (\b -> zip (b_seq b) (SB.unpack $ b_qual b)) brs'
    md' = mk_new_md cigar' md1 seq1 cons_seq

    get_xp br = case M.lookup "XP" (b_exts br) of Just (Int i) -> i ; _ -> 1

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
    [ "F'ing MD field is f'ed up when constructing new MD!"
    , "CIGAR: " ++ show cigs
    , "MD: " ++ show ms
    , "refseq: " ++ show osq
    , "readseq: " ++ show nsq ]


consensus :: Word8 -> [ (Nucleotide, Word8) ] -> (Nucleotide, Word8)
consensus maxq nqs = if qr > 3 then (n0, qr) else (N,0)
  where
    accs :: UArray Nucleotide Int
    accs = accumArray (+) 0 (A,T) [ (n,fromIntegral q) | (n,q) <- nqs, isProperBase n ]

    (n0,q0) : (_,q1) : _ = sortBy (flip $ comparing snd) $ assocs accs
    qr = fromIntegral $ (q0-q1) `min` fromIntegral maxq


