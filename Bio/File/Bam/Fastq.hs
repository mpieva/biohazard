{-# LANGUAGE OverloadedStrings #-}
module Bio.File.Bam.Fastq (
    parseFastq, parseFastq'
                          ) where

-- Parser for FastA/FastQ, @Iteratee@ style, based on Attoparsec.

import Bio.File.Bam
import Bio.Util
import Control.Applicative       hiding ( many )
import Data.Attoparsec.Char8
import Data.Attoparsec.Iteratee
import Data.Bits
import Data.Iteratee             hiding ( length )

import qualified Data.Attoparsec.Char8  as P
import qualified Data.ByteString        as S
import qualified Data.Iteratee.ListLike as I
import qualified Data.Map               as M

-- | Reader for DNA (not protein) sequences in FastA and FastQ.  We read
-- everything vaguely looking like FastA or FastQ, then shoehorn it into
-- a BAM record.  We strive to extract information following
-- more-or-less-established conventions from the header, but we won't
-- support everything under the sun.  Only the canonical variant of
-- FastQ is supported (qualities stored as raw bytes with base 33).
-- Input can be gzipped.
--
-- Supported conventions:
-- * A name suffix of /1 or /2 is turned into the first mate or second
--   mate flag and the read is flagged as paired.
-- * Same for name prefixes of F_ or R_, respectively.
-- * A name prefix of M_ flags the sequence as unpaired and both first
--   and second mate (a merged read).
-- * A name prefix of T_ flags the sequence as unpaired and second mate
--   (an adapter trimmed read).
-- * A name prefix of C_, either before or after any of the other
--   prefixes, is turned into the extra flag XP:i:-1 (result of
--   duplicate removal with unknown depth).
-- * A collection of tags separated from the name by an octothorpe is
--   removed and but into the field ZT (provisional: tags) as text.
--
-- Everything before the first sequence header is ignored.  Headers can
-- start with '>' or '@', we treat both equally.  The first word of the
-- header becomes the read name, a parser for the remainder of the
-- header can be supplied, and the default is to simply ignore it.  The
-- sequence can be split across multiple lines; whitespace, dashes and
-- dots are ignored, IUPAC ambiguity codes are accepted as bases,
-- anything else causes an error.  The sequence ends at a line that is
-- either a header or starts with '+', in the latter case, that line is
-- ignored and must be followed by quality scores.  There must be
-- exactly as many Q-scores as there are bases, followed immediately by
-- a header or end-of-file.  Whitespace is ignored.

parseFastq :: Monad m => Enumeratee S.ByteString [ BamRec ] m a
parseFastq = parseFastq' (const id)

-- | Same as @parseFastq@, but a custom function can be applied to the
-- description string and the parsed record.

parseFastq' :: Monad m => (S.ByteString->BamRec->BamRec) -> Enumeratee S.ByteString [ BamRec ] m a
parseFastq' descr it = do skipJunk ; convStream (parserToIteratee $ (:[]) <$> pRec) it
  where
    isCBase   = inClass "ACGTUBDHVSWMKRYNacgtubdhvswmkryn"
    canSkip c = isSpace c || c == '.' || c == '-'
    isHdr   c = c == '@' || c == '>'

    pRec   = (satisfy isHdr <?> "start marker") *> (makeRecord <$> pName <*> (descr <$> P.takeWhile ('\n' /=)) <*> (pSeq >>= pQual))
    pName  = takeTill isSpace <* skipWhile (\c -> c /= '\n' && isSpace c)  <?> "read name"
    pSeq   =     (:) <$> satisfy isCBase <*> pSeq
             <|> satisfy canSkip *> pSeq 
             <|> pure []                                                   <?> "sequence" 

    pQual sq = (,) sq <$> (char '+' *> skipWhile ('\n' /=) *> pQual' (length sq) <* skipSpace <|> return S.empty)  <?> "qualities"
    pQual' n = S.map (subtract 33) . S.filter (not . isSpace_w8) <$> scan n step
    step 0 _ = Nothing
    step i c | isSpace c = Just i
             | otherwise = Just (i-1)

skipJunk :: Monad m => Iteratee S.ByteString m ()
skipJunk = peek >>= check
  where
    check (Just c) | bad c = I.dropWhile (c2w '\n' /=) >> I.drop 1 >> skipJunk
    check _                = return ()
    bad c = c /= c2w '>' && c /= c2w '@'

makeRecord :: Seqid -> (BamRec->BamRec) -> (String, S.ByteString) -> BamRec
makeRecord name0 extra (sq,qual) = extra $ BamRec {
        b_qname = name,
        b_flag  = flags,
        b_rname = invalidRefseq,
        b_pos   = invalidPos,
        b_mapq  = 0,
        b_cigar = Cigar [],
        b_mrnm  = invalidRefseq,
        b_mpos  = invalidPos,
        b_isize = 0,
        b_seq   = read sq,
        b_qual  = qual,
        b_exts  = tags,
        b_virtual_offset = 0 }
  where
    (name, flags, tags) = checkFR $ checkC $ checkSharp (name0, 0, M.empty)

    checkFR (n,f,t) | "F_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagFirstMate  .|. flagPaired,     t)
                    | "R_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagSecondMate .|. flagPaired,     t)
                    | "M_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagFirstMate  .|. flagSecondMate, t)
                    | "T_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|.                    flagSecondMate, t)
                    | "/1" `S.isSuffixOf` n =        ( rdrop 2 n, f .|. flagFirstMate  .|. flagPaired,     t)
                    | "/2" `S.isSuffixOf` n =        ( rdrop 2 n, f .|. flagSecondMate .|. flagPaired,     t)
                    | otherwise             =        (         n, f,                                       t)

    checkC (n,f,t) | "C_" `S.isPrefixOf` n  = (S.drop 2 n, f, M.insert "XP" (Int (-1)) t)
                   | otherwise              = (         n, f,                          t)

    rdrop n s = S.take (S.length s - n) s

    checkSharp (n,f,t) = case S.split (c2w '#') n of [n',ts] -> (n', f, M.insert "ZT" (Text ts) t)
                                                     _       -> ( n, f,                         t)

----------------------------------------------------------------------------

some_file :: FilePath
some_file = "/mnt/ngs_data/101203_SOLEXA-GA04_00007_PEDi_MM_QF_SR/Ibis/Final_Sequences/s_5_L3280_sequence_merged.txt"

fastq_test :: FilePath -> IO ()
fastq_test = fileDriver $ joinI $ parseFastq $ print_names
            
print_names :: Iteratee [BamRec] IO ()
print_names = I.mapM_ $ S.putStrLn . b_qname

