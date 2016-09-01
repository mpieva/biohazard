module Bio.Align (
    Mode(..),
    myersAlign,
    showAligned
                 ) where

import Bio.Prelude       hiding ( lefts, rights )
import Foreign.C.String         ( CString )
import Foreign.C.Types          ( CInt(..) )
import Foreign.Marshal.Alloc    ( allocaBytes )

import qualified Data.ByteString.Char8      as S
import qualified Data.ByteString.Unsafe     as S
import qualified Data.ByteString.Lazy.Char8 as L

foreign import ccall unsafe "myers_align.h myers_diff" myers_diff ::
        CString -> CInt ->              -- sequence A and length A
        CInt ->                         -- mode (an enum)
        CString -> CInt ->              -- sequence B and length B
        CInt ->                         -- max distance
        CString ->                      -- backtracing space A
        CString ->                      -- backtracing space B
        IO CInt                         -- returns distance

-- | Mode argument for 'myersAlign', determines where free gaps are
-- allowed.
data Mode = Globally  -- ^ align globally, without gaps at either end
          | HasPrefix -- ^ align so that the second sequence is a prefix of the first
          | IsPrefix  -- ^ align so that the first sequence is a prefix of the second
    deriving Enum

-- | Align two strings.  @myersAlign maxd seqA mode seqB@ tries to align
-- @seqA@ to @seqB@, which will work as long as no more than @maxd@ gaps
-- or mismatches are incurred.  The @mode@ argument determines if either
-- of the sequences is allowed to have an overhanging tail.
--
-- The result is the triple of the actual distance (gaps + mismatches)
-- and the two padded sequences.  These sequences are the original
-- sequences with dashes inserted for gaps.
--
-- The algorithm is the O(nd) algorithm by Myers, implemented in C.  A
-- gap and a mismatch score the same.  The strings are supposed to code
-- for DNA, the code understands IUPAC-IUB ambiguity codes.  Two
-- characters match iff there is at least one nucleotide both can code
-- for.  Note that N is a wildcard, while X matches nothing.

myersAlign :: Int -> Bytes -> Mode -> Bytes -> (Int, Bytes, Bytes)
myersAlign maxd seqA mode seqB =
    unsafePerformIO                                 $
    S.unsafeUseAsCStringLen seqA                    $ \(seq_a, len_a) ->
    S.unsafeUseAsCStringLen seqB                    $ \(seq_b, len_b) ->

    -- size of output buffers derives from this:
    -- char *out_a = bt_a + len_a + maxd +2 ;
    -- char *out_b = bt_b + len_b + maxd +2 ;
    allocaBytes (len_a + maxd + 2)                  $ \bt_a ->
    allocaBytes (len_b + maxd + 2)                  $ \bt_b ->

    myers_diff seq_a (fromIntegral len_a)
               (fromIntegral $ fromEnum mode)
               seq_b (fromIntegral len_b)
               (fromIntegral maxd) bt_a bt_b      >>= \dist ->
    if dist < 0
      then return (maxBound, S.empty, S.empty)
      else (,,) (fromIntegral dist) <$>
           S.packCString bt_a <*>
           S.packCString bt_b


-- | Nicely print an alignment.  An alignment is simply a list of
-- strings with inserted gaps to make them align.  We split them into
-- manageable chunks, stack them vertically and add a line showing
-- asterisks in every column where all aligned strings agree.  The
-- result is /almost/ the Clustal format.
showAligned :: Int -> [Bytes] -> [L.ByteString]
showAligned w ss | all S.null ss = []
                 | otherwise = map (L.fromChunks . (:[])) lefts ++
                               L.pack agreement :
                               L.empty :
                               showAligned w rights
  where
    (lefts, rights) = unzip $ map (S.splitAt w) ss
    agreement = map star $ S.transpose lefts
    star str = if S.null str || S.all (== S.head str) str then '*' else ' '

