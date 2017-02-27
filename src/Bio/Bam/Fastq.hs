module Bio.Bam.Fastq ( parseFastq, parseFastq', parseFastqCassava ) where

import Bio.Bam.Header
import Bio.Bam.Rec
import Bio.Prelude hiding ( isSpace )
import Bio.Iteratee
import Data.Attoparsec.ByteString.Char8

import qualified Data.Attoparsec.ByteString.Char8   as P
import qualified Data.ByteString                    as B
import qualified Data.ByteString.Char8              as S
import qualified Data.Vector.Generic                as V

-- ^ Parser for @FastA/FastQ@, 'Iteratee' style, based on
-- "Data.Attoparsec", and written such that it is compatible with module
-- 'Bio.Bam'.  This gives import of @FastA/FastQ@ while respecting some
-- local conventions.

-- | Reader for DNA (not protein) sequences in FastA and FastQ.  We read
-- everything vaguely looking like FastA or FastQ, then shoehorn it into
-- a BAM record.  We strive to extract information following more or
-- less established conventions from the header, but we won't support
-- everything under the sun.  The recognized syntactical warts are
-- converted into appropriate flags and removed.  Only the canonical
-- variant of FastQ is supported (qualities stored as raw bytes with
-- base 33).
--
-- Supported additional conventions:
--
-- * A name suffix of @/1@ or @/2@ is turned into the first mate or second
--   mate flag and the read is flagged as paired.
--
-- * Same for name prefixes of @F_@ or @R_@, respectively.
--
-- * A name prefix of @M_@ flags the sequence as unpaired and merged
--
-- * A name prefix of @T_@ flags the sequence as unpaired and trimmed
--
-- * A name prefix of @C_@, either before or after any of the other
--   prefixes, is turned into the extra flag @XP:i:-1@ (result of
--   duplicate removal with unknown duplicate count).
--
-- * A collection of tags separated from the name by an octothorpe is
--   removed and put into the fields @XI@ and @XJ@ as text.
--
-- * In 'parseFastqCassava' only, if the first word of the description
--   has at least four colon separated subfields, the first if used to
--   flag first/second mate, the second is the \"QC failed\" flag, and
--   the fourth is the index sequence.
--
-- Everything before the first sequence header is ignored.  Headers can
-- start with @\>@ or @\@@, we treat both equally.  The first word of
-- the header becomes the read name, the remainder of the header is
-- ignored.  The sequence can be split across multiple lines;
-- whitespace, dashes and dots are ignored, IUPAC-IUB ambiguity codes
-- are accepted as bases, anything else causes an error.  The sequence
-- ends at a line that is either a header or starts with @\+@, in the
-- latter case, that line is ignored and must be followed by quality
-- scores.  There must be exactly as many Q-scores as there are bases,
-- followed immediately by a header or end-of-file.  Whitespace is
-- ignored.

parseFastq :: Monad m => Enumeratee Bytes [ BamRec ] m a
parseFastq = parseFastq' (const id)

parseFastqCassava :: Monad m => Enumeratee Bytes [ BamRec ] m a
parseFastqCassava = parseFastq' (pdesc . S.split ':' . S.takeWhile (' ' /=))
  where
    pdesc (num:flg:_:idx:_) br = br { b_flag = sum [ if num == "1" then flagFirstMate .|. flagPaired else 0
                                                   , if num == "2" then flagSecondMate .|. flagPaired else 0
                                                   , if flg == "Y" then flagFailsQC else 0
                                                   , b_flag br .&. complement (flagFailsQC .|. flagSecondMate .|. flagPaired) ]
                                    , b_exts = if S.all (`S.elem` "ACGTN") idx then insertE "XI" (Text idx) (b_exts br) else b_exts br }
    pdesc _ br = br

-- | Same as 'parseFastq', but a custom function can be applied to the
-- description string (the part of the header after the sequence name),
-- which can modify the parsed record.  Note that the quality field can
-- end up empty.

{-# WARNING parseFastq' "parseFastq' no longer removes syntactic warts!" #-}
parseFastq' :: Monad m => ( Bytes -> BamRec -> BamRec ) -> Enumeratee Bytes [ BamRec ] m a
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
    pQual' n = B.filter (not . isSpace_w8) <$> scan n step
    step 0 _ = Nothing
    step i c | isSpace c = Just i
             | otherwise = Just (i-1)

skipJunk :: Monad m => Iteratee Bytes m ()
skipJunk = peekStream >>= check
  where
    check (Just c) | bad c = dropWhileStream (c2w '\n' /=) >> dropStream 1 >> skipJunk
    check _                = return ()
    bad c = c /= c2w '>' && c /= c2w '@'

makeRecord :: Seqid -> (BamRec->BamRec) -> (String, Bytes) -> BamRec
makeRecord name extra (sq,qual) = extra $ nullBamRec
        { b_qname = name, b_seq = V.fromList $ read sq, b_qual = V.fromList $ map (Q . subtract 33) $ B.unpack qual }

