module Bio.Bam.Header (
        BamMeta(..),
        parseBamMeta,
        parseBamMetaLine,
        showBamMeta,
        addPG,

        BamKey(..),
        BamHeader(..),
        BamSQ(..),
        BamSorting(..),
        BamOtherShit,

        Refseq(..),
        invalidRefseq,
        isValidRefseq,
        invalidPos,
        isValidPos,
        unknownMapq,
        isKnownMapq,

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
        eflagTrimmed,
        eflagMerged,
        eflagAlternative,
        eflagExactIndex,

        distinctBin,

        MdOp(..),
        readMd,
        showMd
    ) where

import Bio.Prelude           hiding ( uncons )
import Data.ByteString              ( uncons )
import Data.ByteString.Builder
import Data.Sequence                ( (><), (|>) )

import qualified Data.Attoparsec.ByteString.Char8   as P
import qualified Data.ByteString.Char8              as S
import qualified Data.Sequence                      as Z

data BamMeta = BamMeta {
        meta_hdr :: !BamHeader,
        meta_refs :: !Refs,
        meta_other_shit :: [(BamKey, BamOtherShit)],
        meta_comment :: [Bytes]
    } deriving Show

-- | Exactly two characters, for the \"named\" fields in bam.
newtype BamKey = BamKey Word16
    deriving ( Eq, Ord )

instance IsString BamKey where
    {-# INLINE fromString #-}
    fromString [a,b]
        | ord a < 256 && ord b < 256
            = BamKey . fromIntegral $ ord a .|. shiftL (ord b) 8

    fromString s
            = error $ "Not a legal BAM key: " ++ show s

instance Show BamKey where
    show (BamKey a) = [ chr (fromIntegral a .&. 0xff), chr (shiftR (fromIntegral a) 8 .&. 0xff) ]

addPG :: Maybe Version -> IO (BamMeta -> BamMeta)
addPG vn = do
    args <- getArgs
    pn   <- getProgName
    return $ go args pn
  where
    go args pn bm = bm { meta_other_shit = ("PG",pg_line) : meta_other_shit bm }
      where
        pg_line = concat [ [ ("ID", pg_id) ]
                         , [ ("PN", S.pack pn) ]
                         , [ ("CL", S.pack $ unwords args) ]
                         , maybe [] (\v -> [("VN",S.pack (showVersion v))]) vn
                         , map (\p -> ("PP",p)) (take 1 pg_pp)
                         , map (\p -> ("pp",p)) (drop 1 pg_pp) ]

        pg_id : _ = filter (not . flip elem pg_ids) . map S.pack $
                      pn : [ pn ++ '-' : show i | i <- [(1::Int)..] ]

        pg_ids = [ pgid | ("PG",fs) <- meta_other_shit bm, ("ID",pgid) <- fs ]
        pg_pps = [ pgid | ("PG",fs) <- meta_other_shit bm, ("PP",pgid) <- fs ]

        pg_pp  = pg_ids \\ pg_pps


instance Monoid BamMeta where
    mempty = BamMeta mempty noRefs [] []
    a `mappend` b = BamMeta { meta_hdr = meta_hdr a `mappend` meta_hdr b
                            , meta_refs = meta_refs a >< meta_refs b
                            , meta_other_shit = nub $ meta_other_shit a ++ meta_other_shit b
                            , meta_comment = nub $ meta_comment a ++ meta_comment b }

data BamHeader = BamHeader {
        hdr_version :: (Int, Int),
        hdr_sorting :: !BamSorting,
        hdr_other_shit :: BamOtherShit
    } deriving (Show, Eq)

instance Monoid BamHeader where
    mempty = BamHeader (1,0) Unknown []
    a `mappend` b = BamHeader { hdr_version = hdr_version a `min` hdr_version b
                              , hdr_sorting = let u = hdr_sorting a ; v = hdr_sorting b in if u == v then u else Unknown
                              , hdr_other_shit = hdr_other_shit a ++ hdr_other_shit b }

data BamSQ = BamSQ {
        sq_name :: Seqid,
        sq_length :: Int,
        sq_other_shit :: BamOtherShit
    } deriving (Show, Eq)

bad_seq :: BamSQ
bad_seq = BamSQ (error "no SN field") (error "no LN field") []

-- | Possible sorting orders from bam header.  Thanks to samtools, which
-- doesn't declare sorted files properly, we have to have the stupid
-- 'Unknown' state, too.
data BamSorting = Unknown | Unsorted | Grouped | Queryname | Coordinate | GroupSorted
    deriving (Show, Eq)

type BamOtherShit = [(BamKey, Bytes)]

parseBamMeta :: P.Parser BamMeta
parseBamMeta = fixup . foldl' (flip ($)) mempty <$> P.sepBy parseBamMetaLine (P.skipWhile (=='\t') >> P.char '\n')
  where
    fixup meta = meta { meta_other_shit = reverse (meta_other_shit meta)
                      , meta_comment    = reverse (meta_comment    meta) }

parseBamMetaLine :: P.Parser (BamMeta -> BamMeta)
parseBamMetaLine = P.char '@' >> P.choice [hdLine, sqLine, coLine, otherLine]
  where
    hdLine = P.string "HD\t" >>
             (\fns meta -> let fixup hdr = hdr { hdr_other_shit = reverse (hdr_other_shit hdr) }
                           in meta { meta_hdr = fixup $! foldl' (flip ($)) (meta_hdr meta) fns })
               <$> P.sepBy1 (P.choice [hdvn, hdso, hdother]) tabs

    sqLine = P.string "SQ\t" >>
             (\fns meta -> let fixup sq = sq { sq_other_shit = reverse (sq_other_shit sq) }
                               !s = fixup $ foldl' (flip ($)) bad_seq fns
                           in meta { meta_refs = meta_refs meta |> s })
               <$> P.sepBy1 (P.choice [sqnm, sqln, sqother]) tabs

    hdvn = P.string "VN:" >>
           (\a b hdr -> hdr { hdr_version = (a,b) })
             <$> P.decimal <*> ((P.char '.' <|> P.char ':') >> P.decimal)

    hdso = P.string "SO:" >>
           (\s hdr -> hdr { hdr_sorting = s })
             <$> P.choice [ Grouped     <$ P.string "grouped"
                          , Queryname   <$ P.string "queryname"
                          , Coordinate  <$ P.string "coordinate"
                          , GroupSorted <$ P.string "groupsort"
                          , Unsorted    <$ P.string "unsorted"
                          , Unknown     <$ P.skipWhile (\c -> c/='\t' && c/='\n') ]

    sqnm = P.string "SN:" >> (\s sq -> sq { sq_name = s }) <$> pall
    sqln = P.string "LN:" >> (\i sq -> sq { sq_length = i }) <$> P.decimal

    hdother = (\t hdr -> t `seq` hdr { hdr_other_shit = t : hdr_other_shit hdr }) <$> tagother
    sqother = (\t sq  -> t `seq` sq  { sq_other_shit = t : sq_other_shit sq }) <$> tagother

    coLine = P.string "CO\t" >>
             (\s meta -> s `seq` meta { meta_comment = s : meta_comment meta })
               <$> P.takeWhile (/= '\n')

    otherLine = (\k ts meta -> meta { meta_other_shit = (k,ts) : meta_other_shit meta })
                  <$> bamkey <*> (tabs >> P.sepBy1 tagother tabs)

    tagother :: P.Parser (BamKey,Bytes)
    tagother = (,) <$> bamkey <*> (P.char ':' >> pall)

    tabs = P.char '\t' >> P.skipWhile (== '\t')

    pall :: P.Parser Bytes
    pall = P.takeWhile (\c -> c/='\t' && c/='\n')

    bamkey :: P.Parser BamKey
    bamkey = (\a b -> fromString [a,b]) <$> P.anyChar <*> P.anyChar

showBamMeta :: BamMeta -> Builder
showBamMeta (BamMeta h ss os cs) =
    show_bam_meta_hdr h <>
    foldMap show_bam_meta_seq ss <>
    foldMap show_bam_meta_other os <>
    foldMap show_bam_meta_comment cs
  where
    show_bam_meta_hdr (BamHeader (major,minor) so os') =
        byteString "@HD\tVN:" <>
        intDec major <> char7 '.' <> intDec minor <>
        byteString (case so of Unknown     -> mempty
                               Unsorted    -> "\tSO:unsorted"
                               Grouped     -> "\tSO:grouped"
                               Queryname   -> "\tSO:queryname"
                               Coordinate  -> "\tSO:coordinate"
                               GroupSorted -> "\tSO:groupsort") <>
        show_bam_others os'

    show_bam_meta_seq (BamSQ  _  _ []) = mempty
    show_bam_meta_seq (BamSQ nm ln ts) =
        byteString "@SQ\tSN:" <> byteString nm <>
        byteString "\tLN:" <> intDec ln <> show_bam_others ts

    show_bam_meta_comment cm = byteString "@CO\t" <> byteString cm <> char7 '\n'

    show_bam_meta_other (BamKey k,ts) =
        char7 '@' <> word16LE k <> show_bam_others ts

    show_bam_others ts =
        foldMap show_bam_other ts <> char7 '\n'

    show_bam_other (BamKey k,v) =
        char7 '\t' <> word16LE k <> char7 ':' <> byteString v


-- | Reference sequence in Bam
-- Bam enumerates the reference sequences and then sorts by index.  We
-- need to track that index if we want to reproduce the sorting order.
newtype Refseq = Refseq { unRefseq :: Word32 } deriving (Eq, Ord, Ix)

instance Show Refseq where
    showsPrec p (Refseq r) = showsPrec p r

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

-- | The invalid position.
-- Bam uses this value to encode a missing position.
{-# INLINE invalidPos #-}
invalidPos :: Int
invalidPos = -1

-- | Tests whether a position is valid.
-- Returns true unless the the argument equals @invalidPos@.
{-# INLINE isValidPos #-}
isValidPos :: Int -> Bool
isValidPos = (/=) invalidPos

{-# INLINE unknownMapq #-}
unknownMapq :: Int
unknownMapq = 255

isKnownMapq :: Int -> Bool
isKnownMapq = (/=) unknownMapq

-- | A list of reference sequences.
type Refs = Z.Seq BamSQ

-- | The empty list of references.  Needed for BAM files that don't really store alignments.
noRefs :: Refs
noRefs = Z.empty

getRef :: Refs -> Refseq -> BamSQ
getRef refs (Refseq i)
    | i < fromIntegral (Z.length refs) = Z.index refs (fromIntegral i)
    | otherwise                        = BamSQ "*" 0 []


flagPaired, flagProperlyPaired, flagUnmapped, flagMateUnmapped,
 flagReversed, flagMateReversed, flagFirstMate, flagSecondMate,
 flagAuxillary, flagFailsQC, flagDuplicate :: Int

flagPaired         =   0x1
flagProperlyPaired =   0x2
flagUnmapped       =   0x4
flagMateUnmapped   =   0x8
flagReversed       =  0x10
flagMateReversed   =  0x20
flagFirstMate      =  0x40
flagSecondMate     =  0x80
flagAuxillary      = 0x100
flagFailsQC        = 0x200
flagDuplicate      = 0x400

eflagTrimmed, eflagMerged, eflagAlternative, eflagExactIndex :: Int
eflagTrimmed     = 0x1
eflagMerged      = 0x2
eflagAlternative = 0x4
eflagExactIndex  = 0x8


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
compareNames n m = case (uncons n, uncons m) of
        ( Nothing, Nothing ) -> EQ
        ( Just  _, Nothing ) -> GT
        ( Nothing, Just  _ ) -> LT
        ( Just (c,n'), Just (d,m') )
            | is_digit c && is_digit d ->
                let Just (u,n'') = S.readInt n
                    Just (v,m'') = S.readInt m
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


data MdOp = MdNum Int | MdRep Nucleotides | MdDel [Nucleotides] deriving Show

readMd :: Bytes -> Maybe [MdOp]
readMd s | S.null s           = return []
         | isDigit (S.head s) = do (n,t) <- S.readInt s
                                   (MdNum n :) <$> readMd t
         | S.head s == '^'    = let (a,b) = S.break isDigit (S.tail s)
                                in (MdDel (map toNucleotides $ S.unpack a) :) <$> readMd b
         | otherwise          = (MdRep (toNucleotides $ S.head s) :) <$> readMd (S.tail s)

-- | Normalizes a series of 'MdOp's and encodes them in the way BAM and
-- SAM expect it.
showMd :: [MdOp] -> Bytes
showMd = S.pack . flip s1 []
  where
    s1 (MdNum  i : MdNum  j : ms) = s1 (MdNum (i+j) : ms)
    s1 (MdNum  0            : ms) = s1 ms
    s1 (MdNum  i            : ms) = shows i . s1 ms

    s1 (MdRep  r            : ms) = shows r . s1 ms

    s1 (MdDel d1 : MdDel d2 : ms) = s1 (MdDel (d1++d2) : ms)
    s1 (MdDel []            : ms) = s1 ms
    s1 (MdDel ns : MdRep  r : ms) = (:) '^' . shows ns . (:) '0' . shows r . s1 ms
    s1 (MdDel ns            : ms) = (:) '^' . shows ns . s1 ms
    s1 [                        ] = id


-- | Computes the "distinct bin" according to the BAM binning scheme.  If
-- an alignment starts at @pos@ and its CIGAR implies a length of @len@
-- on the reference, then it goes into bin @distinctBin pos len@.
distinctBin :: Int -> Int -> Int
distinctBin beg len = mkbin 14 $ mkbin 17 $ mkbin 20 $ mkbin 23 $ mkbin 26 0
  where end = beg + len - 1
        mkbin n x = if beg `shiftR` n /= end `shiftR` n then x
                    else ((1 `shiftL` (29-n))-1) `div` 7 + (beg `shiftR` n)
