{-# LANGUAGE OverloadedStrings, BangPatterns #-}
module Bio.Bam.Header (
        BamMeta(..),
        parseBamMeta,
        parseBamMetaLine,
        showBamMeta,
        addPG,

        BamHeader(..),
        BamSQ(..),
        BamSorting(..),
        BamOtherShit,

        Refseq(..),
        invalidRefseq,
        isValidRefseq,

        Refs,
        noRefs,
        getRef,

        compareNames,

        flagPaired,
        flagProperlyPaired,
        flagUnmapped,
        flagMateUnmapped,
        flagReversed,
        flagMateReversed,
        flagFirstMate,
        flagSecondMate,
        flagAuxillary,
        flagFailsQC,
        flagDuplicate,
        flagTrimmed,
        flagMerged
    ) where

import Bio.Base
import Control.Applicative
import Data.Ix
import Data.List                    ( (\\) )
import Data.Monoid
import Data.Sequence                ( (<|), (><) )
import Data.Version                 ( Version, showVersion )
import Data.Word                    ( Word32 )
import System.Environment           ( getArgs, getProgName )

import qualified Data.Attoparsec.Char8          as P
import qualified Data.ByteString                as S
import qualified Data.ByteString.Char8          as B
import qualified Data.ByteString.Lazy.Char8     as L
import qualified Data.Foldable                  as F
import qualified Data.Sequence                  as Z

data BamMeta = BamMeta {
        meta_hdr :: BamHeader,
        meta_refs :: Refs,
        meta_other_shit :: [(Char, Char, BamOtherShit)],
        meta_comment :: [S.ByteString]
    } deriving Show


addPG :: Maybe Version -> IO (BamMeta -> BamMeta)
addPG vn = do
    args <- getArgs
    pn   <- getProgName
    return $ go args pn
  where
    go args pn bm = bm { meta_other_shit = ('P','G',pg_line) : meta_other_shit bm }
      where
        pg_line = concat [ [ ('I','D', pg_id) ]
                         , [ ('P','N', B.pack pn) ]
                         , [ ('C','L', B.pack $ unwords args) ]
                         , maybe [] (\v -> [('V','N',B.pack (showVersion v))]) vn
                         , map (\p -> ('P','P',p)) (take 1 pg_pp)
                         , map (\p -> ('p','p',p)) (drop 1 pg_pp) ]

        pg_id : _ = filter (not . flip elem pg_ids) . map B.pack $
                      pn : [ pn ++ '-' : show i | i <- [(1::Int)..] ]

        pg_ids = [ pgid | ('P','G',fs) <- meta_other_shit bm, ('I','D',pgid) <- fs ]
        pg_pps = [ pgid | ('P','G',fs) <- meta_other_shit bm, ('P','P',pgid) <- fs ]

        pg_pp  = pg_ids \\ pg_pps


instance Monoid BamMeta where
    mempty = BamMeta mempty noRefs [] []
    a `mappend` b = BamMeta { meta_hdr = meta_hdr a `mappend` meta_hdr b
                            , meta_refs = meta_refs a >< meta_refs b
                            , meta_other_shit = meta_other_shit a ++ meta_other_shit b
                            , meta_comment = meta_comment a ++ meta_comment b }

data BamHeader = BamHeader {
        hdr_version :: (Int, Int),
        hdr_sorting :: BamSorting,
        hdr_other_shit :: BamOtherShit
    } deriving Show

instance Monoid BamHeader where
    mempty = BamHeader (1,0) Unsorted []
    a `mappend` b = BamHeader { hdr_version = hdr_version a `min` hdr_version b
                              , hdr_sorting = let u = hdr_sorting a ; v = hdr_sorting b in if u == v then u else Unsorted
                              , hdr_other_shit = hdr_other_shit a ++ hdr_other_shit b }

data BamSQ = BamSQ {
        sq_name :: Seqid,
        sq_length :: Int,
        sq_other_shit :: BamOtherShit
    } deriving Show

bad_seq :: BamSQ
bad_seq = BamSQ (error "no SN field") (error "no LN field") []

data BamSorting = Unsorted | Grouped | Queryname | Coordinate | GroupSorted
    deriving (Show, Eq)

type BamOtherShit = [(Char, Char, S.ByteString)]

parseBamMeta :: P.Parser BamMeta
parseBamMeta = foldr ($) mempty <$> P.sepBy parseBamMetaLine (P.skipWhile (=='\t') >> P.char '\n')

parseBamMetaLine :: P.Parser (BamMeta -> BamMeta)
parseBamMetaLine = P.char '@' >> P.choice [hdLine, sqLine, coLine, otherLine]
  where
    hdLine = P.string "HD\t" >>
             (\fns meta -> meta { meta_hdr = foldr ($) (meta_hdr meta) fns })
               <$> P.sepBy1 (P.choice [hdvn, hdso, hdother]) tabs

    sqLine = P.string "SQ\t" >>
             (\fns meta -> meta { meta_refs = foldr ($) bad_seq fns <| meta_refs meta })
               <$> P.sepBy1 (P.choice [sqnm, sqln, sqother]) tabs

    hdvn = P.string "VN:" >>
           (\a b hdr -> hdr { hdr_version = (a,b) })
             <$> P.decimal <*> ((P.char '.' <|> P.char ':') >> P.decimal)

    hdso = P.string "SO:" >>
           (\s hdr -> hdr { hdr_sorting = s })
             <$> P.choice [ Grouped  <$ P.string "grouped"
                          , Queryname <$ P.string "queryname"
                          , Coordinate <$ P.string "coordinate"
                          , GroupSorted <$ P.string "groupsort"
                          , Unsorted <$ P.skipWhile (\c -> c/='\t' && c/='\n') ]

    sqnm = P.string "SN:" >> (\s sq -> sq { sq_name = s }) <$> pall
    sqln = P.string "LN:" >> (\i sq -> sq { sq_length = i }) <$> P.decimal

    hdother = (\t hdr -> hdr { hdr_other_shit = t : hdr_other_shit hdr }) <$> tagother
    sqother = (\t sq  -> sq  { sq_other_shit = t : sq_other_shit sq }) <$> tagother

    coLine = P.string "CO\t" >>
             (\s meta -> meta { meta_comment = s : meta_comment meta })
               <$> P.takeWhile (/= 'n')

    otherLine = (\a b ts meta -> meta { meta_other_shit = (a,b,ts) : meta_other_shit meta })
                  <$> P.anyChar <*> P.anyChar <*> (tabs >> P.sepBy1 tagother tabs)

    tagother :: P.Parser (Char,Char,S.ByteString)
    tagother = (,,) <$> P.anyChar <*> P.anyChar <*> (P.char ':' >> pall)

    tabs = P.char '\t' >> P.skipWhile (== '\t')

    pall :: P.Parser S.ByteString
    pall = P.takeWhile (\c -> c/='\t' && c/='\n')

showBamMeta :: BamMeta -> L.ByteString -> L.ByteString
showBamMeta (BamMeta h ss os cs) =
    show_bam_meta_hdr h .
    F.foldr ((.) . show_bam_meta_seq) id ss .
    foldr ((.) . show_bam_meta_other) id os .
    foldr ((.) . show_bam_meta_comment) id cs
  where
    show_bam_meta_hdr (BamHeader (major,minor) so os') =
        L.append "@HD\tVN:" . L.append (L.pack (show major ++ '.' : show minor)) .
        L.append (case so of Unsorted -> L.empty
                             Grouped  -> "\tSO:grouped"
                             Queryname  -> "\tSO:queryname"
                             Coordinate  -> "\tSO:coordinate"
                             GroupSorted  -> "\tSO:groupsort") .
        show_bam_others os'

    show_bam_meta_seq (BamSQ  _  _ []) = id
    show_bam_meta_seq (BamSQ nm ln ts) =
        L.append "@SQ\tSN:" . L.append (L.fromChunks [nm]) . L.append "\tLN:" .
        L.append (L.pack (show ln)) . show_bam_others ts

    show_bam_meta_comment cm = L.append "@CO\t" . L.append (L.fromChunks [cm]) . L.cons '\n'

    show_bam_meta_other (a,b,ts) =
        L.cons '@' . L.cons a . L.cons b . show_bam_others ts

    show_bam_others ts =
        foldr ((.) . show_bam_other) id ts . L.cons '\n'

    show_bam_other (a,b,v) =
        L.cons '\t' . L.cons a . L.cons b . L.cons ':' . L.append (L.fromChunks [v])

-- | Reference sequence in Bam
-- Bam enumerates the reference sequences and then sorts by index.  We
-- need to track that index if we want to reproduce the sorting order.
newtype Refseq = Refseq { unRefseq :: Word32 } deriving (Show, Eq, Ord, Ix)

instance Enum Refseq where
    succ = Refseq . succ . unRefseq
    pred = Refseq . pred . unRefseq
    toEnum = Refseq . fromIntegral
    fromEnum = fromIntegral . unRefseq
    enumFrom = map Refseq . enumFrom . unRefseq
    enumFromThen (Refseq a) (Refseq b) = map Refseq $ enumFromThen a b
    enumFromTo (Refseq a) (Refseq b) = map Refseq $ enumFromTo a b
    enumFromThenTo (Refseq a) (Refseq b) (Refseq c) = map Refseq $ enumFromThenTo a b c


-- | Tests whether a reference sequence is valid.
-- Returns true unless the the argument equals @invalidRefseq@.
isValidRefseq :: Refseq -> Bool
isValidRefseq = (/=) invalidRefseq

-- | The invalid Refseq.
-- Bam uses this value to encode a missing reference sequence.
invalidRefseq :: Refseq
invalidRefseq = Refseq 0xffffffff

-- | A list of reference sequences.
type Refs = Z.Seq BamSQ

-- | The empty list of references.  Needed for BAM files that don't really store alignments.
noRefs :: Refs
noRefs = Z.empty

getRef :: Refs -> Refseq -> BamSQ
getRef refs (Refseq i)
    | 0 <= i && fromIntegral i <= Z.length refs = Z.index refs (fromIntegral i)
    | otherwise                                 = BamSQ "*" 0 []


flagPaired, flagProperlyPaired, flagUnmapped, flagMateUnmapped, flagReversed, flagMateReversed, flagFirstMate, flagSecondMate,
 flagAuxillary, flagFailsQC, flagDuplicate, flagTrimmed, flagMerged :: Int

flagPaired = 0x1
flagProperlyPaired = 0x2
flagUnmapped = 0x4
flagMateUnmapped = 0x8
flagReversed = 0x10
flagMateReversed = 0x20
flagFirstMate = 0x40
flagSecondMate = 0x80
flagAuxillary = 0x100
flagFailsQC = 0x200
flagDuplicate = 0x400

flagTrimmed = 0x10000
flagMerged  = 0x20000


-- | Compares two sequence names the way samtools does.
-- samtools sorts by "strnum_cmp":
-- . if both strings start with a digit, parse the initial
--   sequence of digits and compare numerically, if equal,
--   continue behind the numbers
-- . else compare the first characters (possibly NUL), if equal
--   continue behind them
-- . else both strings ended and the shorter one counts as
--   smaller (and that part is stupid)

compareNames :: Seqid -> Seqid -> Ordering
compareNames n m = case (S.uncons n, S.uncons m) of
        ( Nothing, Nothing ) -> EQ
        ( Just  _, Nothing ) -> GT
        ( Nothing, Just  _ ) -> LT
        ( Just (c,n'), Just (d,m') )
            | is_digit c && is_digit d ->
                let Just (u,n'') = B.readInt n
                    Just (v,m'') = B.readInt m
                in case u `compare` v of
                    LT -> LT
                    GT -> GT
                    EQ -> n'' `compareNames` m''
            | otherwise -> case c `compare` d of
                    LT -> LT
                    GT -> GT
                    EQ -> n' `compareNames` m'
  where
    is_digit c = 48 <= c && c < 58


