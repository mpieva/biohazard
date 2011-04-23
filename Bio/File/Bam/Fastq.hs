{-# LANGUAGE Rank2Types #-}
module Bio.File.Bam.Fastq where

import Bio.File.Bam
import Control.Applicative
import Data.Char ( toUpper, ord )
import Data.Word ( Word8 )

import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.ByteString.Lazy       as LB
import qualified Data.Map as M

-- | Reader for DNA (not protein) sequences in FastA and FastQ.  We read
-- everything vaguely looking like FastA or FastQ, then shoehorn it into
-- a BAM record.  We strive to extract information following established
-- conventions from the header, but we won't support everything under
-- the sun.  Only the canonical variant of FastQ is supported (qualities
-- stored as raw bytes with base 33).  Input can be gzipped.
--
-- Supported conventions:
-- * A name suffix of /1 or /2 is turned into the first mate or second
--   mate flag and the read is flagged as paired.
-- * Same for name prefixes of F_ or R_.
-- * A name prefix of M_ flags the sequence as unpaired and both first
--   and second mate (a merged read).
-- * A prefix of C_ is turned into the extra flag XD:i:0 (result of
--   duplicate removal with unknown depth).
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

parseFastq :: L.ByteString -> [ BamRec ]
parseFastq = parseFastq' skipDescr
  where skipDescr = many (pTest ('\n' /=)) *> pure id

parseFastq' :: Parser (BamRec->BamRec) -> L.ByteString -> [ BamRec ]
parseFastq' pDescr = run pFasta
  where
    run (P p) = p just1 fail fail
    just1 x y = if L.all isSpace y then x else fail y
    fail s = error $ "parse error near " ++ show (L.take 16 s)

    isHdr c = c == '@' || c == '>'
    isBase c = toUpper c `elem` "ACGTUBDHVSWMKRYN"
    isSpace c = c == '\n' || c == '\r' || c == ' ' || c == '\t'
    canSkip c = isSpace c || c == '.' || c == '-'

    pFasta = many pJunk *> unfold pRec
    pJunk  = pTest (not . isHdr) *> many (pTest ('\n' /=)) *> pTest ('\n' ==)
    pRec :: Parser BamRec
    pRec = pTest isHdr *> (makeRecord <$> pName <*> pDescr <*> (pSeq >>= pQual))
    pName = many (pTest (not . isSpace))
    pSeq = (:) <$> pTest isBase <*> pSeq <|> pTest canSkip *> pSeq <|> pure []
    pQual sq = (,) sq <$> (pSym '+' *> many (pTest ('\n' /=)) *> pQual' (length sq) <|> return [])
    pQual' n = pTest isSpace *> pQual' n <|>
               (if n == 0 then pure [] else cons <$> pGet <*> pQual' (n-1))
                    where cons c qs = (fromIntegral (ord c) - 33) : qs
    unfold p = ((:) <$> p <!> unfold p) <|> pure []

    makeRecord name extra (sq,qual) = extra $ BamRec {
            b_qname = L.pack name,              -- XXX flags?
            b_flag  = 0,                        -- XXX flags?
            b_rname = invalidRefseq,
            b_pos   = invalidPos,
            b_mapq  = 0,
            b_cigar = packCigar [],
            b_mrnm  = invalidRefseq,
            b_mpos  = invalidPos,
            b_isize = 0,
            b_seq   = encodeSeq $ read sq,
            b_qual  = LB.pack qual,
            b_exts  = M.empty,
            b_virtual_offset = 0 }


----------------------------------------------------------------------------
--
-- Module	: HXML.LLParsing
-- Copyright	: (C) 2000-2002 Joe English.  Freely redistributable.
-- License	: "MIT-style"
--
-- Author	: Joe English <jenglish@flightlab.com>
--
-- Simple, non-backtracking, no-lookahead parser combinators.
-- Use with caution!  "Borrowed" and adapted to ByteString input.


{- module LLParsing
    ( pTest , pCheck , pSym ,
    ,(<^>),(<?>)
    , pFoldr , pChainr , pChainl, pTry
    , pRun
    ) where -}

newtype Parser res = P (forall a .
       (res -> L.ByteString -> a)   -- ok continuation
    -> (L.ByteString -> a)          -- failure continuation
    -> (L.ByteString -> a)          -- error continuation
    -> L.ByteString                 -- input
    -> a)                           -- result


pGet :: Parser Char
pGet = P pget
  where
    pget ok f _e i | L.null i = f L.empty
                   | otherwise = ok (L.head i) (L.tail i)

pTest :: (Char -> Bool) -> Parser Char
pTest pred = P (ptest pred)
  where
    ptest p ok f _e l | L.null l     = f L.empty
                      | p (L.head l) = ok (L.head l) (L.tail l)
                      | otherwise    = f l

pSym :: Char -> Parser Char
pSym a = pTest (a==)

pCheck	:: (Char -> Maybe b) -> Parser b
pCheck cmf = P (pcheck cmf)
  where
    pcheck mf ok f _e cs
        | L.null cs = f L.empty
        | otherwise = case mf (L.head cs) of
            Just x	-> ok x (L.tail cs)
            Nothing	-> f cs

pTry :: Parser a -> Parser a
pTry (P pa) = P (\ok f _e i -> pa ok f (\ _i' -> f i) i)

instance Functor Parser where
    fmap f (P pb)	= P (\ok -> pb (ok . f))

instance Applicative Parser where
    pure a = P (\ok _f _e  -> ok a)
    P pa <*> P pb = P (\ok f e -> pa (\a -> pb (ok . a) e e) f e)

instance Alternative Parser where
    empty = P (\_ok f _e -> f)
    (P pa) <|> (P pb) = P (\ok f e -> pa ok (pb ok f e) e)

instance Monad Parser where
    return a = P (\ok _f _e  -> ok a)
    P pa >>= k = P (\ok f e -> pa (\a -> case k a of P pb -> pb ok e e) f e)

infixr 3 <!>
(<!>) :: Parser (a->b) -> Parser a -> Parser b                      -- eager sequence (contains a thinko?)
P pa <!> P pb = P (\ok -> pa (\a i -> ok (a (pb at_eof fail fail i)) L.empty))
  where fail s = error $ "parse error before " ++ show (L.take 16 s)
        at_eof r s = if L.null s then r else fail s

infixl 4 <^>, <?> 

(<?>) :: Parser b -> b -> Parser b                                  -- optional
P pa <?> a = P (\ok _f -> pa ok (ok a))

(<^>) :: Parser b -> Parser c -> Parser (b,c)                       -- sequence
P pa <^> P pb = P (\ok f e -> pa (\a->pb(\b->ok (a,b)) e e) f e)

pFoldr :: (a->b->b) -> b -> Parser a -> Parser b
pFoldr op e p = loop where loop = (op <$> p <*> loop) <?> e

pChainr :: Parser (b -> b -> b) -> Parser b -> Parser b
pChainr op p = loop
  where loop = p <**> ((flip <$> op <*> loop) <?> id)
