{-# LANGUAGE DefaultSignatures, TypeOperators, FlexibleContexts, FlexibleInstances #-}
module Bio.Util.Pretty where

import Bio.Prelude hiding ( Prefix, Infix, Text, (<+>), (<$>) )
import Data.Text.Lazy ( Text )
import GHC.Generics
import Text.PrettyPrint.Leijen.Text hiding ( Pretty(..), (<>), empty )

import qualified Data.Vector.Unboxed as U

-- ^ I can has generic pretty printer?

pprint :: Pretty a => a -> IO ()
pprint = (>> putStrLn []) . displayIO stdout . renderPretty 0.75 132 . pretty 0

prettyList :: Pretty a => [a] -> Doc
prettyList = brackets . align . fillSep . punctuate comma . map (pretty 0)

class Pretty a where
    pretty :: Int -> a -> Doc

    default pretty :: (Generic a, GPretty (Rep a)) => Int -> a -> Doc
    pretty i x = gpretty Pref i (from x)

instance Pretty Double where pretty _ = double

instance (Pretty a, U.Unbox a) => Pretty (U.Vector a) where pretty _ = prettyList . U.toList

data PType = Rec | Pref | Inf Text

appPrec, appPrec1 :: Int
appPrec = 10
appPrec1 = 11

parensIf :: Bool -> Doc -> Doc
parensIf True  = parens
parensIf False = id

class GPretty f where
    gpretty :: PType -> Int -> f a -> Doc

    basic :: f a -> Bool
    basic _ = False

instance GPretty U1 where
    gpretty _ _ _ = mempty
    basic _ = True

instance Pretty c => GPretty (K1 i c) where
    gpretty _ i (K1 fp) = pretty i fp

instance (GPretty f, Constructor c) => GPretty (M1 C c f) where
    gpretty t i c@(M1 fp)
        | conIsRecord c = text (fromString (conName c)) <+>
                          braces (align (gpretty Rec i fp))
        | otherwise = case conFixity c of
            Prefix -> parensIf (i > appPrec && not (basic fp)) $
                      text (fromString (conName c)) <+>
                       gpretty t appPrec1 fp
            Infix _ m -> parensIf (i > m) $
                         braces (align (gpretty t m fp))

instance (GPretty f, Selector c) => GPretty (M1 S c f) where
    gpretty t i c@(M1 fp) = case t of
        Pref  -> gpretty t i fp
        Inf _ -> gpretty t i fp
        Rec   -> text (fromString (selName c)) <+>
                 char '=' <+> align (gpretty t i fp)

instance GPretty f => GPretty (M1 D c f) where
    gpretty t i (M1 fp) = gpretty t i fp

instance (GPretty f, GPretty g) => GPretty (f :+: g) where
    gpretty t i (L1 fp) = gpretty t i fp
    gpretty t i (R1 fp) = gpretty t i fp

instance (GPretty f, GPretty g) => GPretty (f :*: g) where
    gpretty t@Rec     n (a :*: b) = gpretty t   n   a <> char ',' <$> gpretty t   n b
    gpretty t@Pref    n (a :*: b) = gpretty t (n+1) a             </> gpretty t (n+1) b
    gpretty t@(Inf s) n (a :*: b) = gpretty t   n   a <> text   s </> gpretty t   n b

