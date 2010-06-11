module Bio.File.TwoBit (
    TwoBitFile,
    openTwoBit,
    closeTwoBit,
    withTwoBit,

    getSubseq
) where

{-

Would you believe it?  The 2bit format stores blocks of Ns in a table at
the beginning of a sequence, then packs four bases into a byte.  So it
is neither possible nor necessary to store Ns in the main sequence, and
you would think they aren't stored there, right?  And they aren't.
Instead Ts are stored which the reader has to replace with Ns.  

How stupid is that?

-}

import Bio.Base

import           Control.Exception
import           Control.Monad
import           Data.Array.Unboxed
import           Data.Bits
import           Data.Binary.Get
import qualified Data.ByteString.Char8 as S
import qualified Data.ByteString.Lazy as L
import qualified Data.ByteString.Lazy.Char8 as LC
import qualified Data.IntMap as I
import           Data.IORef
import qualified Data.Map as M
import           Data.Maybe
import           Data.Word
import           Numeric
import           System.IO
import           System.IO.Unsafe

data TwoBitFile = TBF {
    tbf_handle :: !Handle,
    tbf_get_word32 :: !(Get Int),
    tbf_seqs :: !(M.Map Seqid (IORef TwoBitSequence))
}

data TwoBitSequence = Untouched { tbs_offset :: !Int }
                    | Indexed   { tbs_s_blocks :: !( M.Map Int Int )
                                , tbs_m_blocks :: !( I.IntMap Int )
                                , tbs_dna_offset :: !Int }

openTwoBit :: FilePath -> IO TwoBitFile
openTwoBit fp = do
    h <- openFile fp ReadMode
    raw <- pGetContents h 0
    ~(g,ix) <- return . flip runGet raw $ do
                sig <- getWord32be
                getWord32 <- case sig of
                        0x1A412743 -> return $ fromIntegral `fmap` getWord32be
                        0x4327411A -> return $ fromIntegral `fmap` getWord32le
                        _          -> fail $ "invalid signature " ++ showHex sig []

                
                version <- getWord32
                unless (version == 0) $ fail $ "wrong version " ++ show version

                nseqs <- getWord32
                getWord32

                (,) getWord32 `fmap` replicateM nseqs ( liftM2 (,)
                        ( getWord8 >>= getByteString . fromIntegral )
                        ( getWord32 >>= return . Untouched ) )

    m <- let foldM' [    ] acc = return acc
             foldM' (x:xs) acc = cons x acc >>= foldM' xs
             cons (k,v) m = v `seq` newIORef v >>= \r -> return $! M.insert (S.copy k) r m
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
                                                -- mb <- read_block_list
                                                skip_block_list
                                                len <- getWord32 >> bytesRead

                                                return $! Indexed (M.fromList $ to_good_blocks ds nb)
                                                                  (I.empty {-I.fromList mb-}) (ofs + fromIntegral len)
                                   writeIORef r $! sq'
                                   return sq'
  where
    getWord32 = tbf_get_word32 tbf
    skip_block_list = getWord32 >>= skip . (*) 8
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

do_frag :: Int -> Int -> Sense -> M.Map Int Int -> (Integer -> IO L.ByteString) -> Int -> IO [Nucleotide]
do_frag start len revcomplp s_blocks raw ofs0 = do
    dna <- get_dna (case revcomplp of Forward -> fwd_nt ; Reverse -> cmp_nt) 
                   start len final_blocks raw ofs0
    return $ case revcomplp of { Forward -> dna ; Reverse -> reverse dna }

  where
    (left_junk, mfirst, left_clipped) = M.splitLookup start s_blocks

    left_fragment | M.null left_junk = Nothing
                  | otherwise  = case M.findMax left_junk of
        (start1, len1) -> 
            let d = start - start1 
            in if d >= len1 then Nothing
                            else Just (start, len1-d)

    (right_clipped, _) = M.split (start+len) .
                         maybe id (uncurry M.insert) left_fragment .
                         maybe id (M.insert start) mfirst $ left_clipped

    right_fragment | M.null right_clipped = Nothing
                   | otherwise = case M.findMax right_clipped of
        (startn, lenn) ->
            let l' = start + len - startn 
            in if l' <= 0 then Nothing
                          else Just (startn, min l' lenn)

    rdrop1 [] = [] ; rdrop1 [_] = [] ; rdrop1 (x:xs) = x : rdrop1 xs

    final_blocks = rdrop1 (M.toAscList right_clipped) ++ maybe [] (:[]) right_fragment

    fwd_nt = (!!) [T,C,A,G] . fromIntegral
    cmp_nt = (!!) [A,G,T,C] . fromIntegral



get_dna :: (Word8 -> Nucleotide) -> Int -> Int -> [(Int, Int)] -> (Integer -> IO L.ByteString) -> Int -> IO [Nucleotide]
get_dna _nt _start total [] _raw _ofs0              = return $ replicate total N
get_dna _nt _start total _  _raw _ofs0 | total <= 0 = return []
get_dna  nt  start total blocks0@((start1, len1) : blocks) raw ofs
    | start /= start1 =
        (replicate (start1-start) N ++) `fmap`
        get_dna nt start1 (total-start1+start) blocks0 raw ofs
    | otherwise = do
        s <- L.take bytes' `fmap` raw (fromIntegral $ ofs + (ofs' `div` 4))
        let dna = take len1 . drop (start1 .&. 3) . L.foldr toDNA [] $ s
        (++) dna `fmap` get_dna nt (start+len1) (total-len1) blocks raw ofs
  where ofs'   = start1 .&. complement 3
        bytes' = fromIntegral $ (len1 + (start1 .&. 3) + 3) `div` 4
        toDNA b = (++) [ nt (3 .&. (b `shiftR` x)) | x <- [6,4,2,0] ]



getSubseq :: TwoBitFile -> Range -> IO [Nucleotide]
getSubseq tbf (Range { r_pos = Pos { p_seq = chr, p_start = start, p_sense = strand }, r_length = len }) = do
             ref <- maybe (fail $ S.unpack chr ++ " doesn't exist") return 
                    $ M.lookup chr (tbf_seqs tbf)
             sq1 <- read_block_index tbf ref
             do_frag start len strand (tbs_s_blocks sq1) (pGetContents $ tbf_handle tbf) (tbs_dna_offset sq1)

test :: IO ()
test =
    withTwoBit "/mnt/sequencedb/mammalian_genomes/9606.2bit"
            (\f -> getSubseq f (Range (Pos (S.pack "chr1") Reverse 217280) 60000) >>= print)

pGetContents :: Handle -> Integer -> IO L.ByteString
pGetContents hdl ofs = L.fromChunks `fmap` go ofs
  where
    chunk_size = 4096
    go o = unsafeInterleaveIO $ liftM2 (:)
            (hSeek hdl AbsoluteSeek o >> S.hGet hdl chunk_size)
            (go $ o + fromIntegral chunk_size)

