{-# LANGUAGE BangPatterns #-}
-- Re-wrap alignments to obey the given length of the reference
-- sequence.
--
-- The idea is that a circular reference sequence has been extended
-- artificially to facilitate alignment.  Now the declared length in the
-- header is wrong, and the alignments overhang the end.  Here we split
-- those alignments into two, one for the beginning, one for the end of
-- the sequence, then soft-mask out the inappropriate parts.
--
-- What's the best course of action, operationally?  As usual, we need
-- to decide whether to rely on sorted input and whether to produce
-- sorted output, and how much to copy senselessly.
--
-- In a sane world, this program runs precisely once, after alignment,
-- and output is piped somewhere.  So that's what we do:  input is
-- unsorted, so is output, output is piped (and hence uncompressed).
-- We also fix the header while we're at it.
--
-- We try to fix the map quality for the affected reads as follows:  if
-- a read has map quality 0 (meaning multiple equally good hits), we
-- check the XA field.  If it reports exactly one additional alignment,
-- and it matches the primary alignment when transformed to canonical
-- coordinates, we remove XA and set MAPQ to 37.

import Bio.Bam
import Bio.Base
import Control.Monad                    ( when )
import Data.Foldable                    ( toList )
import Data.Version                     ( showVersion )
import Paths_biohazard                  ( version )
import System.Environment               ( getArgs, getProgName )
import System.Exit                      ( exitFailure )
import System.IO                        ( hPutStr )

import qualified Data.ByteString.Char8  as S
import qualified Data.Map               as M
import qualified Data.Sequence          as Z

usage :: IO a
usage = do pn <- getProgName
           hPutStr stderr $ pn ++ ", version " ++ showVersion version ++
                "\nUsage: " ++ pn ++ " [chrom:length...]\n\
                \Pipes a BAM file from stdin to stdout and for every 'chrom'\n\
                \mentioned on the command line, wraps alignments to a new \n\
                \target length of 'length'.\n"
           exitFailure

main :: IO ()
main = getArgs >>= \args ->
       when (null args) usage >>= \_ ->
       enumHandle defaultBufSize stdin >=> run $
       joinI $ decodeAnyBam $ \hdr -> do
           add_pg <- liftIO (addPG $ Just version)
           let (ltab, seqs') = parseArgs (meta_refs hdr) args
           joinI $ mapChunks (concatMap (rewrap (M.fromList ltab)))
                 $ protectTerm $ pipeRawBamOutput (add_pg hdr { meta_refs = seqs' })

parseArgs :: Refs -> [String] -> ([(Refseq,(Int,S.ByteString))], Refs)
parseArgs refs | Z.null refs = error $ "no target sequences found (empty input?)"
               | otherwise   = foldl parseArg ([],refs)
  where
    parseArg (sqs, h) arg = case break (==':') arg of
        (nm,':':r) -> case reads r of
            [(l,[])] | l > 0 -> case filter (S.isPrefixOf (S.pack nm) . sq_name . snd) $ zip [0..] $ toList h of
                [(k,a)] | sq_length a >= l -> ( (Refseq $ fromIntegral k,(l, sq_name a)):sqs, Z.update k (a { sq_length = l }) h )
                        | otherwise -> error $ "cannot wrap " ++ show nm ++ " to " ++ show l
                                            ++ ", which is more than the original " ++ show (sq_length a)
                [] -> error $ "no match for target sequence " ++ show nm
                _ -> error $ "target sequence " ++ show nm ++ " is ambiguous"
            _ -> error $ "couldn't parse length " ++ show r ++ " for " ++ show nm
        _ -> error $ "couldn't parse argument " ++ show arg



-- | This runs both stages of the rewrapping: First normalize alignments
-- (POS must be in the canonical interval) and fix XA, MPOS, MAPQ where
-- appropriate, then duplicate the read and softmask the noncanonical
-- parts.  Rmdup fits in between the two, hence the split
rewrap :: M.Map Refseq (Int,S.ByteString) -> BamRaw -> [BamRaw]
rewrap m br = maybe [br] (\(l,nm) -> wrapTo l $ normalizeTo nm l br)
              $ M.lookup (br_rname br) m
