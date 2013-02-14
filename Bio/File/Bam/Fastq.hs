-- | Parser for FastA/FastQ, @Iteratee@ style, based on Attoparsec, and
-- written such that it is compatible with module @Bio.File.Bam@.

{-# LANGUAGE OverloadedStrings #-}
module Bio.File.Bam.Fastq (
    parseFastq, parseFastq', removeWarts
                          ) where

import Bio.File.Bam
import Bio.Iteratee
import Control.Applicative       hiding ( many )
import Data.Attoparsec.Char8
import Data.Attoparsec.Iteratee
import Data.Bits

import qualified Data.Attoparsec.Char8  as P
import qualified Data.ByteString        as B
import qualified Data.ByteString.Char8  as S
import qualified Data.Iteratee.ListLike as I
import qualified Data.Map               as M

-- | Reader for DNA (not protein) sequences in FastA and FastQ.  We read
-- everything vaguely looking like FastA or FastQ, then shoehorn it into
-- a BAM record.  We strive to extract information following
-- more-or-less-established conventions from the header, but we won't
-- support everything under the sun.  The recognized syntactical warts
-- are converted into appropriate flags and removed.  Only the canonical
-- variant of FastQ is supported (qualities stored as raw bytes with
-- base 33).  Input can be gzipped.
--
-- Supported conventions:
-- * A name suffix of /1 or /2 is turned into the first mate or second
--   mate flag and the read is flagged as paired.
-- * Same for name prefixes of F_ or R_, respectively.
-- * A name prefix of M_ flags the sequence as unpaired and merged
-- * A name prefix of T_ flags the sequence as unpaired and trimmed
-- * A name prefix of C_, either before or after any of the other
--   prefixes, is turned into the extra flag XP:i:-1 (result of
--   duplicate removal with unknown depth).
-- * A collection of tags separated from the name by an octothorpe is
--   removed and put into the fields XI and YI as text.
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
    pQual' n = B.map (subtract 33) . B.filter (not . isSpace_w8) <$> scan n step
    step 0 _ = Nothing
    step i c | isSpace c = Just i
             | otherwise = Just (i-1)

skipJunk :: Monad m => Iteratee S.ByteString m ()
skipJunk = I.peek >>= check
  where
    check (Just c) | bad c = I.dropWhile (c2w '\n' /=) >> I.drop 1 >> skipJunk
    check _                = return ()
    bad c = c /= c2w '>' && c /= c2w '@'

makeRecord :: Seqid -> (BamRec->BamRec) -> (String, S.ByteString) -> BamRec
makeRecord name extra (sq,qual) = extra $ removeWarts $ nullBamRec
        { b_qname = name, b_seq = read sq, b_qual = qual }

-- | Remove syntactic warts from old read names.
removeWarts :: BamRec -> BamRec
removeWarts br = br { b_qname = name, b_flag = flags, b_exts = tags }
  where
    (name, flags, tags) = checkFR $ checkC $ checkSharp (b_qname br, b_flag br, b_exts br)

    checkFR (n,f,t) | "F_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagFirstMate  .|. flagPaired, t)
                    | "R_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagSecondMate .|. flagPaired, t)
                    | "M_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagMerged,                    t)
                    | "T_" `S.isPrefixOf` n = checkC (S.drop 2 n, f .|. flagTrimmed,                   t)
                    | "/1" `S.isSuffixOf` n =        ( rdrop 2 n, f .|. flagFirstMate  .|. flagPaired, t)
                    | "/2" `S.isSuffixOf` n =        ( rdrop 2 n, f .|. flagSecondMate .|. flagPaired, t)
                    | otherwise             =        (         n, f,                                   t)

    checkC (n,f,t) | "C_" `S.isPrefixOf` n  = (S.drop 2 n, f, M.insert "XP" (Int (-1)) t)
                   | otherwise              = (         n, f,                          t)

    rdrop n s = S.take (S.length s - n) s

    checkSharp (n,f,t) = case S.split '#' n of [n',ts] -> (n', f, insertTags ts t)
                                               _       -> ( n, f,               t)

    insertTags ts t | S.null y  = M.insert "XI" (Text ts) t
                    | otherwise = M.insert "XI" (Text  x) $ M.insert "YI" (Text $ S.tail y) t
        where (x,y) = S.break (== ',') ts

----------------------------------------------------------------------------

some_file :: FilePath
some_file = "/mnt/ngs_data/101203_SOLEXA-GA04_00007_PEDi_MM_QF_SR/Ibis/Final_Sequences/s_5_L3280_sequence_merged.txt"

fastq_test :: FilePath -> IO ()
fastq_test = fileDriver $ joinI $ parseFastq $ print_names
            
print_names :: Iteratee [BamRec] IO ()
print_names = I.mapM_ $ S.putStrLn . b_qname

