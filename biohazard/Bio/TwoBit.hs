module Bio.TwoBit (
        module Bio.Base,

        TwoBitFile,
        openTwoBit,
        closeTwoBit,
        withTwoBit,

        getSubseq,
        getSubseqAscii,
        getSeqnames,
        hasSequence,
        getSeqLength,
        clampPosition,
        getRandomSeq,

        Mask(..)
    ) where

{-
Would you believe it?  The 2bit format stores blocks of Ns in a table at
the beginning of a sequence, then packs four bases into a byte.  So it
is neither possible nor necessary to store Ns in the main sequence, and
you would think they aren't stored there, right?  And they aren't.
Instead Ts are stored which the reader has to replace with Ns.

How stupid is that?

The sensible way to treat these is probably to just say there are two
kinds of implied annotation (repeats and large gaps for a typical
genome), which can be interpreted in whatever way fits.

TODO:  use Judy for the Int->Int mappings?  Or sorted arrays with binary
       search?
TODO:  use 'Iteratee's for the IO instead of lazy reading
-}

import           Bio.Base
import           Control.Applicative
import           Control.Exception
import           Control.Monad
import           Data.Bits
import           Data.Binary.Get
import qualified Data.ByteString.Char8 as S
import qualified Data.ByteString.Lazy as L
import           Data.Char (toLower)
import qualified Data.IntMap as I
import           Data.IORef
import qualified Data.Map as M
import           Data.Maybe
import           Numeric
import           System.IO
import           System.IO.Unsafe
import           System.Random

data TwoBitFile = TBF {
    tbf_handle :: !Handle,
    tbf_get_word32 :: !(Get Int),
    tbf_seqs :: !(M.Map Seqid (IORef TwoBitSequence))
}

data TwoBitSequence = Untouched { _tbs_offset    :: {-# UNPACK #-} !Int }
                    | Indexed   { tbs_n_blocks   :: !(I.IntMap Int)
                                , tbs_m_blocks   :: !(I.IntMap Int)
                                , tbs_dna_offset :: {-# UNPACK #-} !Int
                                , tbs_dna_size   :: {-# UNPACK #-} !Int }

openTwoBit :: FilePath -> IO TwoBitFile
openTwoBit fp = do
    h <- openFile fp ReadMode
    raw <- pGetContents h 0
    let (g,ix) = flip runGet raw $ do
                sig <- getWord32be
                getWord32 <- case sig of
                        0x1A412743 -> return $ fromIntegral `fmap` getWord32be
                        0x4327411A -> return $ fromIntegral `fmap` getWord32le
                        _          -> fail $ "invalid .2bit signature " ++ showHex sig []


                version <- getWord32
                unless (version == 0) $ fail $ "wrong .2bit version " ++ show version

                nseqs <- getWord32
                _reserved <- getWord32

                (,) getWord32 `fmap` repM nseqs ( liftM2 (,)
                        ( getWord8 >>= getLazyByteString . fromIntegral )
                        ( liftM Untouched getWord32 ) )

    TBF h g <$> loop ix M.empty
  where
    loop [        ] m = return m
    loop ((k,v):xs) m = do r <- newIORef $! v ; loop xs $! M.insert (shelve k) r m



closeTwoBit :: TwoBitFile -> IO ()
closeTwoBit = hClose . tbf_handle

withTwoBit :: FilePath -> (TwoBitFile -> IO a) -> IO a
withTwoBit f = bracket (openTwoBit f) closeTwoBit

readBlockIndex :: TwoBitFile -> IORef TwoBitSequence -> IO TwoBitSequence
readBlockIndex tbf r = do
    sq <- readIORef r
    case sq of Indexed {} -> return sq
               Untouched ofs -> do c <- pGetContents (tbf_handle tbf) (fromIntegral ofs)
                                   let sq' = flip runGet c $ do
                                                ds <- getWord32
                                                nb <- readBlockList
                                                mb <- readBlockList
                                                len <- getWord32 >> bytesRead

                                                return $! Indexed (I.fromList nb) (I.fromList mb)
                                                                  (ofs + fromIntegral len) ds
                                   writeIORef r $! sq'
                                   return sq'
  where
    getWord32 = tbf_get_word32 tbf
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

getFwdSubseqWith :: (Integer -> IO L.ByteString) -> Int         -- reader fn, dna offset
                 -> I.IntMap Int -> I.IntMap Int                -- N blocks, M blocks
                 -> (Word8 -> Mask -> a)                        -- mask function
                 -> Int -> Int                                  -- start, len
                 -> IO [a]                                      -- result
getFwdSubseqWith raw ofs n_blocks m_blocks nt start len =
    do_mask (takeOverlap start n_blocks `mergeblocks` takeOverlap start m_blocks) start .
    take len . drop (start .&. 3) . L.foldr toDNA [] <$>
    raw (fromIntegral $ ofs + (start `shiftR` 2))
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
getSubseqWith :: (Nucleotide -> Mask -> a) -> TwoBitFile -> Range -> IO [a]
getSubseqWith maskf tbf (Range { r_pos = Pos { p_seq = chr, p_start = start }, r_length = len }) = do
    ref <- maybe (fail $ S.unpack chr ++ " doesn't exist") return $ M.lookup chr (tbf_seqs tbf)
    sq1 <- readBlockIndex tbf ref

    let go = getFwdSubseqWith (pGetContents $ tbf_handle tbf) (tbs_dna_offset sq1)
                              (tbs_n_blocks sq1) (tbs_m_blocks sq1)

    if start < 0 then reverse <$> go (maskf . cmp_nt) (-start-len) len
                 else             go (maskf . fwd_nt)  start      len
  where
    fwd_nt = (!!) [nucT, nucC, nucA, nucG] . fromIntegral
    cmp_nt = (!!) [nucA, nucG, nucT, nucC] . fromIntegral


-- | Extract a subsequence with typical masking:  soft masking is
-- ignored, hard masked regions are replaced with Ns.
getSubseq :: TwoBitFile -> Range -> IO [Nucleotide]
getSubseq = getSubseqWith mymask
  where
    mymask n None = n
    mymask n Soft = n
    mymask _ Hard = nucN
    mymask _ Both = nucN

-- | Extract a subsequence with masking for biologists:  soft masking is
-- done by lowercasing, hard masking by printing an N.
getSubseqAscii :: TwoBitFile -> Range -> IO String
getSubseqAscii = getSubseqWith mymask
  where
    mymask n None = showNucleotide n
    mymask n Soft = toLower (showNucleotide n)
    mymask _ Hard = 'N'
    mymask _ Both = 'N'


pGetContents :: Handle -> Integer -> IO L.ByteString
pGetContents hdl ofs = L.fromChunks `fmap` go ofs
  where
    chunk_size = 32000
    go o = unsafeInterleaveIO $ liftM2 (:)
            (hSeek hdl AbsoluteSeek o >> S.hGet hdl chunk_size)
            (go $ o + fromIntegral chunk_size)

getSeqnames :: TwoBitFile -> [Seqid]
getSeqnames = M.keys . tbf_seqs

hasSequence :: TwoBitFile -> Seqid -> Bool
hasSequence tbf sq = isJust . M.lookup sq . tbf_seqs $ tbf

getSeqLength :: TwoBitFile -> Seqid -> IO Int
getSeqLength tbf chr = do
             ref <- maybe (fail $ shows chr " doesn't exist") return
                    $ M.lookup chr (tbf_seqs tbf)
             sq1 <- readBlockIndex tbf ref
             return $ tbs_dna_size sq1

-- | limits a range to a position within the actual sequence
clampPosition :: TwoBitFile -> Range -> IO Range
clampPosition g (Range (Pos n start) len) = do
    size <- getSeqLength g n

    let start' = if start < 0 then max start (-size) else start
        end'   = min (start + len) $ if start < 0 then 0 else size
    return $ Range (Pos n start') (end' - start')


getRandomSeq :: TwoBitFile -> IO (([Nucleotide] -> Bool) -> Int -> IO (Range, [Nucleotide]))
getRandomSeq tbf = do
    let names = getSeqnames tbf
    lengths <- mapM (getSeqLength tbf) names
    let total = sum lengths
    let frags = I.fromList $ zip (scanl (+) 0 lengths) names

    let draw good l = do p <- randomRIO (1,total)
                         d <- randomRIO (False,True)
                         let Just ((o,s),_) = I.maxViewWithKey $ fst $ I.split p frags
                         r' <- clampPosition tbf $ Range (Pos s (p-o)) l
                         sq <- getSubseq tbf $ if d then r' else reverseRange r'
                         if r_length r' == l && good sq
                           then return (r', sq)
                           else draw good l
    return draw

