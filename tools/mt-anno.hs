import Anno
import Seqs
import Xlate

import Bio.Align
import Control.Applicative
import Data.Char
import Prelude
import Text.Printf

import qualified Data.ByteString.Char8 as S

main :: IO ()
main = do
    smps <- fromFasta . S.lines . S.filter (/= '\r') <$> S.getContents

    let (tabs, fas) = unzip $ map (uncurry do_anno) smps

    putStrLn $ unlines $ concatMap (++ [[]]) tabs
    putStrLn $ unlines $ concatMap (++ [[]]) fas


do_anno :: S.ByteString -> S.ByteString -> ([String], [String])
do_anno smp_name raw_sample = (tab, fa)
  where
    (_, rCRS, sample) = myersAlign 3000 {- maxd -} raw_rCRS Globally raw_sample
    xpose1 = xpose rCRS sample
    xpose_anno g = g { start = xpose1 (start g), end = xpose1 (end g) }

    tab = to_tab (S.unpack smp_name) $ map xpose_anno $ rCRS_anno
    fa  = concatMap (to_fasta smp_name sample xpose1) rCRS_anno


to_fasta :: S.ByteString -> S.ByteString -> (Int -> Int) -> Anno -> [String]
to_fasta smp_name smp f Gene{..} = case what of CDS -> go ; CDS' -> go ; _ -> []
  where
    go = let s' = f start
             e' = f end
             prot = case init $ get_protein smp (s',e') of
                        'I' : rest -> 'M' : rest ; x -> x
             hdr = printf ">%s [gene=%s] [protein=%s] [location=%s]"
                          (S.unpack smp_name) name prod loc
             loc | s' <= e'  = shows s' ".." ++ show e'
                 | otherwise = "complement(" ++ shows e' ".." ++ shows s' ")"
         in hdr : chunk prot

    chunk s = case splitAt 70 s of (u,v) | null v -> [u]
                                         | otherwise -> u : chunk v


fromFasta :: [S.ByteString] -> [(S.ByteString, S.ByteString)]
fromFasta ls = case dropWhile (not . isHeader) ls of
    [      ] -> []
    (h:rest) -> case break isHeader rest of
        (body,rest') -> (S.drop 1 h, S.map toUpper $ S.concat body) : fromFasta rest'
  where
    isHeader s = not (S.null s) && S.head s == '>'

