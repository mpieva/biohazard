{-# LANGUAGE BangPatterns #-}
{-# OPTIONS_GHC -Wall #-}
module Xlate where

import Data.List

import qualified Data.Map as M
import qualified Data.ByteString.Char8 as S
import qualified Data.IntMap as I

import Seqs

-- aligned sequences in, coodinate on first in, coordinate on second out
xpose :: S.ByteString -> S.ByteString -> Int -> Int
xpose ref smp = \p -> I.findWithDefault (-1) p m
  where
    (!m,_,_) = foldl' advance (I.empty, 0, 0) $ S.zip ref smp
    advance (!m,!p1,!p2) (r,s) = let !p1' = if r == '-' then p1 else 1+p1
                                     !p2' = if s == '-' then p2 else 1+p2
                                 in if r == '-' then (m,p1',p2')
                                    else (I.insert p1' p2' m, p1', p2')

-- diffz :: CDS -> [(String, Int, Char, Char)]
-- diffz cds@(CDS _ nm _) = [ (nm, i, r, b) | (i,r,b) <- zip3 [1..] aa_ref aa_bnt, r /= b ]
  -- where (aa_ref, aa_bnt) = get_protein cds

get_protein :: S.ByteString -> (Int,Int) -> String
get_protein ns (s,e) = translate $ cutout
  where
    cutout | s <= e = (take (e-s+1) $ drop (s-1) $ filter (/= '-') $ S.unpack ns) ++ "AA"
           | otherwise = (map compl $ reverse $
                          take (s-e+1) $ drop (e-1) $ filter (/= '-') $ S.unpack ns) ++ "AA"

    compl 'A' = 'T'
    compl 'C' = 'G'
    compl 'G' = 'C'
    compl 'T' = 'A'


translate :: String -> String
translate (a:b:c:s) = m : translate s
    where m = M.findWithDefault (error $ show (a,b,c)) (a,b,c) mito_code
translate _ = []

standard_code :: M.Map (Char,Char,Char)  Char
standard_code = M.fromList $ zip3 base1 base2 base3 `zip` aas
  where
    aas   = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

mito_code :: M.Map (Char,Char,Char)  Char
mito_code = M.insert ('A','G','A') '*' $
            M.insert ('A','G','G') '*' $
            M.insert ('A','T','A') 'M' $
            M.insert ('T','G','A') 'W' $ standard_code

