{-# LANGUAGE TypeOperators, ScopedTypeVariables #-}
module Bio.Util.Pretty where

import Bio.Prelude                       hiding ( Prefix, Infix, (<+>), (<$>) )
import Data.Text.Encoding                       ( decodeUtf8With )
import Data.Text.Encoding.Error                 ( lenientDecode )
import Data.Text.Lazy                           ( fromStrict )
import GHC.Generics
import Text.PrettyPrint.Leijen.Text             ( (<+>), (<$>), (</>) )

import qualified Data.Attoparsec.Text           as A
import qualified Data.Vector                    as V
import qualified Data.Vector.Unboxed            as U
import qualified Text.PrettyPrint.Leijen.Text   as P

-- ^ I can has generic pretty printer?

pprint :: Pretty a => a -> IO ()
pprint = (>> putStrLn []) . P.displayIO stdout . P.renderPretty 0.75 132 . pretty 0

pshow :: Pretty a => a -> LazyText
pshow = P.displayT . P.renderPretty 0.75 132 . pretty 0

pparse :: Parse a => Text -> Either String a
pparse = A.parseOnly (parse 0)

prettyList :: Pretty a => [a] -> P.Doc
prettyList = P.brackets . P.align . P.fillSep . P.punctuate P.comma . map (pretty 0)

parseList :: Parse a => A.Parser [a]
parseList = a_brackets $ parse 0 `A.sepBy` a_comma

a_brackets :: A.Parser a -> A.Parser a
a_brackets p = A.char '[' *> A.skipSpace *> p <* A.skipSpace <* A.char ']' <* A.skipSpace

a_parens :: A.Parser a -> A.Parser a
a_parens p = A.char '(' *> A.skipSpace *> p <* A.skipSpace <* A.char ')' <* A.skipSpace

a_comma :: A.Parser Char
a_comma = A.char ',' <* A.skipSpace

class Pretty a where
    pretty :: Int -> a -> P.Doc

default_pretty :: (Generic a, GPretty (Rep a)) => Int -> a -> P.Doc
default_pretty i x = gpretty Pref i (from x)

class Parse a where
    parse :: Int -> A.Parser a

default_parse :: (Generic a, GParse (Rep a)) => Int -> A.Parser a
default_parse i = to `fmap` gparse Pref i

instance Pretty    Int where pretty _ = P.int
instance Pretty Double where pretty _ = P.double
instance Pretty  Bytes where pretty _ = P.text . fromStrict . decodeUtf8With lenientDecode

instance Parse    Int where parse _ = A.signed A.decimal
instance Parse Double where parse _ = A.double

instance (Pretty a, Pretty b) => Pretty (a,b) where
    pretty p (a,b) = P.parens $ pretty p a <> P.char ',' <+> pretty p b

instance (Parse  a, Parse  b) => Parse  (a,b) where
    parse  p = a_parens $ pure (,) <*> parse p <* A.char ',' <* A.skipSpace <*> parse p

instance Pretty a => Pretty [a] where pretty _ = prettyList
instance Parse  a => Parse  [a] where parse  _ = parseList

instance Pretty a => Pretty (V.Vector a) where pretty _ = prettyList . V.toList
instance Parse a  => Parse  (V.Vector a) where parse  _ = V.fromList `fmap` parseList

instance (Pretty a, U.Unbox a) => Pretty (U.Vector a) where pretty _ = prettyList . U.toList
instance (Parse a,  U.Unbox a) => Parse  (U.Vector a) where parse  _ = U.fromList `fmap` parseList

data PType = Rec | Pref | Inf String

appPrec, appPrec1 :: Int
appPrec = 10
appPrec1 = 11

parensIf :: Bool -> P.Doc -> P.Doc
parensIf True  = P.parens
parensIf False = id

class GPretty f where
    gpretty :: PType -> Int -> f a -> P.Doc

    basic :: f a -> Bool
    basic _ = False

class GParse f where
    gparse :: PType -> Int -> A.Parser (f a)

    basic' :: f a -> Bool
    basic' _ = False

instance GPretty U1 where
    gpretty _ _ _ = mempty
    basic _ = True

instance GParse U1 where
    gparse _ _ = return U1
    basic' _ = True


instance Pretty c => GPretty (K1 i c) where
    gpretty _ i (K1 fp) = pretty i fp

instance Parse c => GParse (K1 i c) where
    gparse _ i = K1 `fmap` parse i


instance (GPretty f, Constructor c) => GPretty (M1 C c f) where
    gpretty t i c@(M1 fp)
        | conIsRecord c = parensIf (i > appPrec && not (basic fp)) $
                          P.text (fromString (conName c)) <+>
                          P.braces (P.align (gpretty Rec i fp))
        | otherwise = case conFixity c of
            Prefix -> parensIf (i > appPrec && not (basic fp)) $
                      P.text (fromString (conName c)) <+>
                      gpretty t appPrec1 fp
            Infix _ m -> parensIf (i > m) $ P.align (gpretty t m fp)

-- XXX What to do with the precedence thingy?
-- XXX Both the prefix and infix forms should be tried.
instance (GParse f, Constructor c) => GParse (M1 C c f) where
    gparse t i | conIsRecord c = fmap M1 $ a_parensIf (i > appPrec1 && not (basic' fp)) $
                                           record_form <|> prefix_form
               | otherwise = fmap M1 $ case conFixity c of
                    Prefix    -> a_parensIf (i > appPrec1 && not (basic' fp)) $ prefix_form
                    Infix _ _ -> a_parensIf (i > appPrec1 && not (basic' fp)) $ prefix_form
      where
        record_form = A.string (fromString (conName c)) *> A.skipSpace *> a_braces (gparse Rec i)
        prefix_form = A.string (fromString (conName c)) *> A.skipSpace *> gparse t appPrec1

        c = M1 undefined :: M1 C c f a
        fp = undefined :: f a

-- XXX What to do with the precedence thingy?
a_parensIf :: Bool -> A.Parser a -> A.Parser a
a_parensIf b p = A.char '(' *> A.skipSpace *> p <* A.skipSpace <* A.char ')' <* A.skipSpace
                 <|> if b then empty else p

a_braces :: A.Parser a -> A.Parser a
a_braces p = A.char '{' *> A.skipSpace *> p <* A.skipSpace <* A.char '}' <* A.skipSpace

instance (GPretty f, Selector c) => GPretty (M1 S c f) where
    gpretty t i c@(M1 fp) = case t of
        Pref  -> gpretty t i fp
        Inf _ -> gpretty t i fp
        Rec   -> P.text (fromString (selName c)) <+>
                 P.char '=' <+> P.align (gpretty t i fp)

instance (GParse f, Selector c) => GParse (M1 S c f) where
    gparse t i = fmap M1 $ gparse t i <|> A.string (fromString (selName c)) *> A.skipSpace
                                          *> A.char '=' *> A.skipSpace *> gparse t i
      where
        c = M1 undefined :: M1 S c f a

instance GPretty f => GPretty (M1 D c f) where
    gpretty t i (M1 fp) = gpretty t i fp

instance GParse f => GParse (M1 D c f) where
    gparse t i = fmap M1 (gparse t i)

instance (GPretty f, GPretty g) => GPretty (f :+: g) where
    gpretty t i (L1 fp) = gpretty t i fp
    gpretty t i (R1 fp) = gpretty t i fp

instance (GParse f, GParse g) => GParse (f :+: g) where
    gparse t i = fmap L1 (gparse t i) <|> fmap R1 (gparse t i)

instance (GPretty f, GPretty g) => GPretty (f :*: g) where
    gpretty t@Rec     n (a :*: b) = gpretty t   n   a <> P.char ','            <$> gpretty t   n b
    gpretty t@Pref    n (a :*: b) = gpretty t (n+1) a                          </> gpretty t (n+1) b
    gpretty t@(Inf s) n (a :*: b) = gpretty t   n   a <> P.text (fromString s) </> gpretty t   n b

instance (GParse f, GParse g) => GParse (f :*: g) where
    gparse t@Rec     n = pure (:*:) <*> gparse t n
                                    <*  A.skipSpace
                                    <*  a_comma
                                    <*> gparse t n
    gparse t@Pref    n = pure (:*:) <*> gparse t (n+1)
                                    <*  A.skipSpace
                                    <*> gparse t (n+1)
    gparse t@(Inf s) n = pure (:*:) <*> gparse t n
                                    <*  A.skipSpace
                                    <*  A.string (fromString s)
                                    <*  A.skipSpace
                                    <*> gparse t n

