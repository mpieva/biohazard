{-# LANGUAGE BangPatterns #-}
module Bio.TwoBit (
        module Bio.Base,

        TwoBitFile,
        openTwoBit,

        getSubseq,
        getSubseqWith,
        getSubseqAscii,
        getSubseqMasked,
        getSeqnames,
        hasSequence,
        getSeqLength,
        clampPosition,
        getRandomSeq,

        Mask(..)
    ) where

import           Bio.Base
import           Control.Applicative
import           Control.Monad
import           Data.Bits
import           Data.Binary.Get
import qualified Data.ByteString as B
import qualified Data.ByteString.Lazy as L
import           Data.Char (toLower)
import qualified Data.IntMap as I
import qualified Data.Map as M
import           Data.Maybe
import           Numeric
import           System.IO.Posix.MMap
import           System.Random

-- ^ Would you believe it?  The 2bit format stores blocks of Ns in a table at
-- the beginning of a sequence, then packs four bases into a byte.  So it
-- is neither possible nor necessary to store Ns in the main sequence, and
-- you would think they aren't stored there, right?  And they aren't.
-- Instead Ts are stored which the reader has to replace with Ns.
--
-- The sensible way to treat these is probably to just say there are two
-- kinds of implied annotation (repeats and large gaps for a typical
-- genome), which can be interpreted in whatever way fits.  And that's why
-- we have 'Mask' and 'getSubseqWith'.
--
-- TODO:  use binary search for the Int->Int mappings?

data TwoBitFile = TBF {
    tbf_raw :: B.ByteString,
    tbf_seqs :: !(M.Map Seqid TwoBitSequence)
}

data TwoBitSequence = Indexed { tbs_n_blocks   :: !(I.IntMap Int)
                              , tbs_m_blocks   :: !(I.IntMap Int)
                              , tbs_dna_offset :: {-# UNPACK #-} !Int
                              , tbs_dna_size   :: {-# UNPACK #-} !Int }

-- | Brings a 2bit file into memory.  The file is mmap'ed, so it will
-- not work on streams that are not actual files.  It's also unsafe if
-- the file is modified in any way.
openTwoBit :: FilePath -> IO TwoBitFile
openTwoBit fp = do raw <- unsafeMMapFile fp
                   return $ flip runGet (L.fromChunks [raw]) $ do
                            sig <- getWord32be
                            getWord32 <- case sig of
                                    0x1A412743 -> return $ fromIntegral `fmap` getWord32be
                                    0x4327411A -> return $ fromIntegral `fmap` getWord32le
                                    _          -> fail $ "invalid .2bit signature " ++ showHex sig []

                            version <- getWord32
                            unless (version == 0) $ fail $ "wrong .2bit version " ++ show version

                            nseqs <- getWord32
                            _reserved <- getWord32

                            TBF raw <$> foldM (\ix _ -> do !key <- getWord8 >>= getByteString . fromIntegral
                                                           !off <- getWord32
                                                           return $! M.insert key (mkBlockIndex raw getWord32 off) ix
                                              ) M.empty [1..nseqs]


mkBlockIndex :: B.ByteString -> Get Int -> Int -> TwoBitSequence
mkBlockIndex raw getWord32 ofs = runGet getBlock $ L.fromChunks [B.drop ofs raw]
  where
    getBlock = do ds <- getWord32
                  nb <- readBlockList
                  mb <- readBlockList
                  len <- getWord32 >> bytesRead
                  return $! Indexed (I.fromList nb) (I.fromList mb) (ofs + fromIntegral len) ds

    readBlockList = getWord32 >>= \n -> liftM2 zip (repM n getWord32) (repM n getWord32)

-- | Repeat monadic action 'n' times.  Returns result in reverse(!) order.
repM :: Monad m => Int -> m a -> m [a]
repM n0 m = go [] n0
  where
    go acc 0 = return acc
    go acc n = m >>= \x -> x `seq` go (x:acc) (n-1)

takeOverlap :: Int -> I.IntMap Int -> [(Int,Int)]
takeOverlap k m = dropWhile far_left $
                  maybe id (\(kv,_) -> (:) kv) (I.maxViewWithKey left) $
                  maybe id (\v -> (:) (k,v)) middle $
                  I.toAscList right
  where
    (left, middle, right) = I.splitLookup k m
    far_left (s,l) = s+l <= k

data Mask = None | Soft | Hard | Both deriving (Eq, Ord, Enum, Show)

getFwdSubseqWith :: B.ByteString -> Int                         -- raw data, dna offset
                 -> I.IntMap Int -> I.IntMap Int                -- N blocks, M blocks
                 -> (Word8 -> Mask -> a)                        -- mask function
                 -> Int -> Int -> [a]                           -- start, len, result
getFwdSubseqWith raw ofs n_blocks m_blocks nt start len =
    do_mask (takeOverlap start n_blocks `mergeblocks` takeOverlap start m_blocks) start .
    take len . drop (start .&. 3) .
    B.foldr toDNA [] .
    B.take (len `div` 4 +1) .
    B.drop (fromIntegral $ ofs + (start `shiftR` 2)) $ raw
  where
    toDNA b = (++) [ 3 .&. (b `shiftR` x) | x <- [6,4,2,0] ]

    do_mask            _ _ [] = []
    do_mask [          ] _ ws = map (`nt` None) ws
    do_mask ((s,l,m):is) p ws = let l0 = (p-s) `max` 0 in
                                map (`nt` None) (take l0 ws) ++
                                map (`nt`    m) (take l (drop l0 ws)) ++
                                do_mask is (p+l0+l) (drop (l+l0) ws)

-- | Merge blocks of Ns and blocks of Ms into single list of blocks with
-- masking annotation.  Gaps remain.  Used internally only.
mergeblocks :: [(Int,Int)] -> [(Int,Int)] -> [(Int,Int,Mask)]
mergeblocks ((_,0):nbs) mbs = mergeblocks nbs mbs
mergeblocks nbs ((_,0):mbs) = mergeblocks nbs mbs

mergeblocks ((ns,nl):nbs) ((ms,ml):mbs)
    | ns < ms   = let l = min (ms-ns) nl in (ns,l, Hard) : mergeblocks ((ns+l,nl-l):nbs) ((ms,ml):mbs)
    | ms < ns   = let l = min (ns-ms) ml in (ms,l, Soft) : mergeblocks ((ns,nl):nbs) ((ms+l,ml-l):mbs)
    | otherwise = let l = min nl ml in (ns,l, Both) : mergeblocks ((ns+l,nl-l):nbs) ((ms+l,ml-l):mbs)

mergeblocks ((ns,nl):nbs) [] = (ns,nl, Hard) : mergeblocks nbs []
mergeblocks [] ((ms,ml):mbs) = (ms,ml, Soft) : mergeblocks [] mbs

mergeblocks [     ] [     ] = []


-- | Extract a subsequence and apply masking.  TwoBit file can represent
-- two kinds of masking (hard and soft), where hard masking is usually
-- realized by replacing everything by Ns and soft masking is done by
-- lowercasing.  Here, we take a user supplied function to apply
-- masking.
getSubseqWith :: (Nucleotide -> Mask -> a) -> TwoBitFile -> Range -> [a]
getSubseqWith maskf tbf (Range { r_pos = Pos { p_seq = chr, p_start = start }, r_length = len }) = do
    let sq1 = maybe (error $ unpackSeqid chr ++ " doesn't exist") id $ M.lookup chr (tbf_seqs tbf)
    let go = getFwdSubseqWith (tbf_raw tbf) (tbs_dna_offset sq1) (tbs_n_blocks sq1) (tbs_m_blocks sq1)
    if start < 0
        then reverse $ go (maskf . cmp_nt) (-start-len) len
        else           go (maskf . fwd_nt)   start      len
  where
    fwd_nt = (!!) [nucT, nucC, nucA, nucG] . fromIntegral
    cmp_nt = (!!) [nucA, nucG, nucT, nucC] . fromIntegral


-- | Extract a subsequence without masking.
getSubseq :: TwoBitFile -> Range -> [Nucleotide]
getSubseq = getSubseqWith (\n _ -> n)

-- | Extract a subsequence with typical masking:  soft masking is
-- ignored, hard masked regions are replaced with Ns.
getSubseqMasked :: TwoBitFile -> Range -> [Nucleotides]
getSubseqMasked = getSubseqWith mymask
  where
    mymask (N n) None = Ns $ 1 `shiftL` fromIntegral n
    mymask (N n) Soft = Ns $ 1 `shiftL` fromIntegral n
    mymask    _  Hard = nucsN
    mymask    _  Both = nucsN

-- | Extract a subsequence with masking for biologists:  soft masking is
-- done by lowercasing, hard masking by printing an N.
getSubseqAscii :: TwoBitFile -> Range -> String
getSubseqAscii = getSubseqWith mymask
  where
    mymask n None = showNucleotide n
    mymask n Soft = toLower (showNucleotide n)
    mymask _ Hard = 'N'
    mymask _ Both = 'N'


getSeqnames :: TwoBitFile -> [Seqid]
getSeqnames = M.keys . tbf_seqs

hasSequence :: TwoBitFile -> Seqid -> Bool
hasSequence tbf sq = isJust . M.lookup sq . tbf_seqs $ tbf

getSeqLength :: TwoBitFile -> Seqid -> Int
getSeqLength tbf chr =
    maybe (error $ shows chr " doesn't exist") tbs_dna_size $
    M.lookup chr (tbf_seqs tbf)

-- | limits a range to a position within the actual sequence
clampPosition :: TwoBitFile -> Range -> Range
clampPosition tbf (Range (Pos n start) len) = Range (Pos n start') (end' - start')
  where
    size   = getSeqLength tbf n
    start' = if start < 0 then max start (-size) else start
    end'   = min (start + len) $ if start < 0 then 0 else size


-- | Sample a piece of random sequence uniformly from the genome.
-- Only pieces that are not hard masked are sampled, soft masking is
-- allowed, but not reported.
-- On a 32bit platform, this will fail for genomes larger than 1G bases.
-- However, if you're running this code on a 32bit platform, you have
-- bigger problems to worry about.
getRandomSeq :: RandomGen g => TwoBitFile                   -- ^ 2bit file
                            -> Int                          -- ^ desired length
                            -> g                            -- ^ RNG
                            -> ((Range, [Nucleotide]), g)   -- ^ position, sequence, new RNG
getRandomSeq tbf len = draw
  where
    names = getSeqnames tbf
    lengths = map (getSeqLength tbf) names
    total = sum lengths
    frags = I.fromList $ zip (scanl (+) 0 lengths) names

    draw g0 | good      = ((r', sq), gn)
            | otherwise = draw gn
      where
        (p0, gn) = randomR (0, 2*total-1) g0
        p = p0 `shiftR` 1
        Just ((o,s),_) = I.maxViewWithKey $ fst $ I.split (p+1) frags
        r' = clampPosition tbf $ Range (Pos s (p-o)) len
        sq = catMaybes $ getSubseqWith mask2maybe tbf $ if odd p0 then r' else reverseRange r'
        good = r_length r' == len && length sq == len

        mask2maybe n None = Just n
        mask2maybe n Soft = Just n
        mask2maybe n Hard = Nothing
        mask2maybe n Both = Nothing

