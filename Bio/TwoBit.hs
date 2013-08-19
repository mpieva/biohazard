module Bio.TwoBit (
        module Bio.Base,

        TwoBitFile,
        openTwoBit,
        closeTwoBit,
        withTwoBit,

        getSubseq,
        getSeqnames,
        hasSequence,
        getSeqLength,
        clampPosition,
        getRandomSeq
    ) where

-- TODO: proper masking is unsupported right now

{-
Would you believe it?  The 2bit format stores blocks of Ns in a table at
the beginning of a sequence, then packs four bases into a byte.  So it
is neither possible nor necessary to store Ns in the main sequence, and
you would think they aren't stored there, right?  And they aren't.
Instead Ts are stored which the reader has to replace with Ns.

How stupid is that?

The sensible way to treat these is probably to just say there are two
kinds of implied annotation (repeats and large gaps for a typical
genome), which can be interpreted in whatever way fits.  All of this
isn't really supported right now.

Note to self:  use Judy for the Int->Int mappings?  Or (gasp!) sorted
arrays with binary search?
-}

import           Bio.Base
import           Control.Exception
import           Control.Monad
import           Data.Array.Unboxed
import           Data.Bits
import           Data.Binary.Get
import qualified Data.ByteString.Char8 as S
import qualified Data.ByteString.Lazy as L
import qualified Data.IntMap as I
import           Data.IORef
import qualified Data.Map as M
import           Data.Maybe
import           Data.Word
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
                    | Indexed   { tbs_s_blocks   :: I.IntMap Int
                                , tbs_m_blocks   :: I.IntMap Int
                                , tbs_dna_offset :: {-# UNPACK #-} !Int
                                , tbs_dna_size   :: {-# UNPACK #-} !Int }

openTwoBit :: FilePath -> IO TwoBitFile
openTwoBit fp = do
    h <- openFile fp ReadMode
    raw <- pGetContents h 0
    ~(g,ix) <- return . flip runGet raw $ do
                sig <- getWord32be
                getWord32 <- case sig of
                        0x1A412743 -> return $ fromIntegral `fmap` getWord32be
                        0x4327411A -> return $ fromIntegral `fmap` getWord32le
                        _          -> fail $ "invalid .2bit signature " ++ showHex sig []


                version <- getWord32
                unless (version == 0) $ fail $ "wrong .2bit version " ++ show version

                nseqs <- getWord32
                _reserved <- getWord32

                (,) getWord32 `fmap` replicateM nseqs ( liftM2 (,)
                        ( getWord8 >>= getLazyByteString . fromIntegral )
                        ( getWord32 >>= return . Untouched ) )

    m <- let foldM' [    ] acc = return acc
             foldM' (x:xs) acc = cons x acc >>= foldM' xs
             cons (k,v) m = v `seq` newIORef v >>= \r -> return $! M.insert (shelve k) r m
         in foldM' ix M.empty
    return $! TBF h g m

closeTwoBit :: TwoBitFile -> IO ()
closeTwoBit = hClose . tbf_handle

withTwoBit :: FilePath -> (TwoBitFile -> IO a) -> IO a
withTwoBit f = bracket (openTwoBit f) closeTwoBit

read_block_index :: TwoBitFile -> IORef TwoBitSequence -> IO TwoBitSequence
read_block_index tbf r = do
    sq <- readIORef r
    case sq of Indexed {} -> return sq
               Untouched ofs -> do c <- pGetContents (tbf_handle tbf) (fromIntegral ofs)
                                   let sq' = flip runGet c $ do
                                                ds <- getWord32
                                                nb <- read_block_list
                                                mb <- read_block_list
                                                len <- getWord32 >> bytesRead

                                                return $! Indexed (I.fromList $ to_good_blocks ds nb)
                                                                  (I.fromList mb)
                                                                  (ofs + fromIntegral len) ds
                                   writeIORef r $! sq'
                                   return sq'
  where
    getWord32 = tbf_get_word32 tbf
    read_block_list = getWord32 >>= \n -> liftM2 zip (read_word_list n) (read_word_list n)
    read_word_list n = listArray (0,n-1) `fmap` repM n getWord32 >>= \arr ->
                       (arr :: UArray Int Int) `seq` return (elems arr)

to_good_blocks :: Int -> [(Int,Int)] -> [(Int,Int)]
to_good_blocks total = tgb 0
  where
    tgb p [              ] = ( p, total - p ) : []
    tgb p ((start,l):rest) = ( p, start - p ) : tgb (start+l) rest


repM :: Monad m => Int -> m a -> m [a]
repM 0 _ = return []
repM n m = m >>= \x -> seq x (repM (n-1) m >>= return . (x:))

do_frag :: Int -> Int -> Bool -> I.IntMap Int -> (Integer -> IO L.ByteString) -> Int -> IO [Nucleotide]
do_frag start0 len revcomplp s_blocks raw ofs0 = do
    dna <- get_dna (if revcomplp then cmp_nt else fwd_nt)
                   start len final_blocks raw ofs0
    return $ if revcomplp then reverse dna else dna

  where
    start = if revcomplp then start0 - len else start0
    (left_junk, mfirst, left_clipped) = I.splitLookup start s_blocks

    left_fragment = case I.maxViewWithKey left_junk of
        Nothing -> Nothing
        Just ((start1, len1), _) ->
            let d = start - start1
            in if d >= len1 then Nothing
                            else Just (start, len1-d)

    (right_clipped, _) = I.split (start+len) .
                         maybe id (uncurry I.insert) left_fragment .
                         maybe id (I.insert start) mfirst $ left_clipped

    right_fragment = case I.maxViewWithKey right_clipped of
        Nothing -> Nothing
        Just ((startn, lenn), _) ->
            let l' = start + len - startn
            in if l' <= 0 then Nothing
                          else Just (startn, min l' lenn)

    rdrop1 [] = [] ; rdrop1 [_] = [] ; rdrop1 (x:xs) = x : rdrop1 xs

    final_blocks = rdrop1 (I.toAscList right_clipped) ++ maybe [] (:[]) right_fragment

    fwd_nt = (!!) [nucT, nucC, nucA, nucG] . fromIntegral
    cmp_nt = (!!) [nucA, nucG, nucT, nucC] . fromIntegral



get_dna :: (Word8 -> Nucleotide) -> Int -> Int -> [(Int, Int)] -> (Integer -> IO L.ByteString) -> Int -> IO [Nucleotide]
get_dna _nt _start total [] _raw _ofs0              = return $ replicate total nucN
get_dna _nt _start total _  _raw _ofs0 | total <= 0 = return []
get_dna  nt  start total blocks0@((start1, len1) : blocks) raw ofs
    | start /= start1 =
        (replicate (start1-start) nucN ++) `fmap`
        get_dna nt start1 (total-start1+start) blocks0 raw ofs
    | otherwise = do
        s <- L.take bytes' `fmap` raw (fromIntegral $ ofs + (ofs' `div` 4))
        let dna = take len1 . drop (start1 .&. 3) . L.foldr toDNA [] $ s
        (++) dna `fmap` get_dna nt (start+len1) (total-len1) blocks raw ofs
  where ofs'   = start1 .&. complement 3
        bytes' = fromIntegral $ (len1 + (start1 .&. 3) + 3) `div` 4
        toDNA b = (++) [ nt (3 .&. (b `shiftR` x)) | x <- [6,4,2,0] ]



getSubseq :: TwoBitFile -> Range -> IO [Nucleotide]
getSubseq tbf (Range { r_pos = Pos { p_seq = chr, p_start = start }, r_length = len }) = do
             ref <- maybe (fail $ S.unpack chr ++ " doesn't exist") return
                    $ M.lookup chr (tbf_seqs tbf)
             sq1 <- read_block_index tbf ref
             let go | start < 0 = do_frag (-start-len) len True
                    | otherwise = do_frag   start      len False
             go (tbs_s_blocks sq1) (pGetContents $ tbf_handle tbf) (tbs_dna_offset sq1)

pGetContents :: Handle -> Integer -> IO L.ByteString
pGetContents hdl ofs = L.fromChunks `fmap` go ofs
  where
    chunk_size = 4096
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
             sq1 <- read_block_index tbf ref
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

