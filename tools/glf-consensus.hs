{-# LANGUAGE BangPatterns #-}
import Control.Applicative ( (<$>) )
import Control.Monad
import Control.Monad.CatchIO
import Control.Monad.Trans.Class
import Data.Char ( isSpace, toLower, chr )
import Data.List ( intercalate, sort )
import System.Console.GetOpt
import System.IO
import System.Environment ( getArgs, getProgName )
import System.Exit

import qualified Data.ByteString.Char8 as S
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.Map as M

import qualified Data.Iteratee.ListLike as I

import Bio.Base
import Bio.Glf
import Bio.Iteratee

data Config = Config {
    conf_min_qual :: Int,
    conf_call     :: [Int] -> [(Int, Char)],
    conf_output   :: Iteratee String IO (),
    conf_input    :: GlfInput,
    conf_conv     :: Formatter,
    conf_mkname   :: S.ByteString -> String }

type GlfInput = (GlfSeq -> Enumeratee [GlfRec] String IO ())
             -> (S.ByteString -> Enumerator String IO ())
             -> Enumerator String IO ()

options :: [ OptDescr (Config -> IO Config) ]
options = [
    Option "1" ["haploid"]
        (NoArg (\c -> return $ c { conf_call = haploid_call }))
        "Force haploid consensus",
    Option "2" ["diploid"]
        (NoArg (\c -> return $ c { conf_call = diploid_call }))
        "Allow diploid consensus",
    Option "m" ["min-qual"]
        (ReqArg (\a c -> readIO a >>= \m -> return $ c { conf_min_qual = m }) "Q")
        "Require minimum quality of Q",
    Option "o" ["output"]
        (ReqArg (\fp c -> return $ c { conf_output = iterToFile fp }) "FILE")
        "Write output to FILE instead of stdout",
    Option "q" ["fastq"]
        (NoArg (\c -> return $ c { conf_conv = print_fastq }))
        "Write FastQ instead of FastA",
    Option "I" ["identifier"]
        (ReqArg (\n c -> return $ c { conf_mkname = subst_name n }) "ID")
        "Use ID as identifier for consensus",
    Option "if" ["input"]
        (ReqArg (\fp c -> return $ c { conf_input = enum_glf_file fp }) "FILE")
        "Read input from FILE instead of stdin",
    Option "h?" ["help", "usage"]
        (NoArg (usage exitSuccess))
        "Print this help" ]


usage :: IO a -> Config -> IO a
usage e _ = getProgName >>= \p -> putStrLn (usageInfo (blurb p) options) >> e
  where blurb prg =
            "Usage: " ++ prg ++ " [Option...] [FastA-File...]\n" ++
            "Reads GLF from stdin and prints the contained consensus sequence in\n" ++
            "FastA/FastQ format.  Gaps are filled with a reference sequence if known\n" ++
            "from the FastA files on the command line, otherwise with Ns."

iterToFile :: FilePath -> Iteratee String IO ()
iterToFile fp = bracket (lift $ openFile fp WriteMode)
                        (lift . hClose)
                        (mapChunksM_ . hPutStr)

defaultConfig :: Config
defaultConfig = Config 0 diploid_call (mapChunksM_ putStr) (enum_glf_handle stdin) print_fasta S.unpack

main :: IO ()
main = do (opts, files, errors) <- getOpt Permute options <$> getArgs
          unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
          Config min_qual call output input conv mkname <- foldl (>>=) (return defaultConfig) opts
          refs <- M.fromList . concatMap readFasta <$> mapM L.readFile files

          hPutStrLn stderr $
                "known reference sequences: [" ++ intercalate ", "
                [ show (L.unpack k) ++ " (" ++ show (L.length v) ++ ")" | (k,v) <- M.toList refs ]
                ++ "]"

          let per_file :: Seqid -> Enumerator String IO ()
              per_file _genome_name = return

              per_seq :: GlfSeq -> Enumeratee [GlfRec] String IO ()
              per_seq glfseq = extract1consensus (mkRef refs glfseq) call min_qual
                               ><> conv (mkname $ glf_seqname glfseq)

          input per_seq per_file output >>= run

-- get the "most likely" consensus, defined as:
-- - as many reference bases or else Ns as were skipped from the previous record, then
-- - if there's an insert, the most likely insert sequence (may be empty)
-- - if there's a deletion, skip the most likely number of bases (may be zero)
-- - else the most likely base

mkRef :: M.Map L.ByteString L.ByteString -> GlfSeq -> Int -> Int -> QSeq
mkRef refs glfseq = case M.lookup (L.fromChunks [glf_seqname glfseq]) refs of
                Nothing -> \o l -> replicate (min l (glf_seqlen glfseq - o)) ('N',2)
                Just s  -> \o l -> let l' = fromIntegral $ min l (glf_seqlen glfseq - o)
                                   in [ (toLower b,30) | b <- L.unpack $ L.take l' $ L.drop (fromIntegral o) s ]

type QSeq = [(Char,Int)]    -- sequence w/ quality

extract1consensus :: Monad m
                  => (Int -> Int -> QSeq)
                  -> ([Int] -> [(Int,Char)])           -- call function
                  -> Int                               -- minimum quality
                  -> Enumeratee [GlfRec] QSeq m r      -- eats records, emits calls
extract1consensus ref call min_qual oit = liftI $ scan oit 0 0
  where
    -- rec_pos: position of last record
    -- ref_pos: first position in reference we haven't handled
    scan k        !_ !ref_pos (EOF        x) = lift  $ enumPure1Chunk (ref ref_pos maxBound) >=> enumChunk (EOF x) $ k
    scan k !rec_pos_ !ref_pos (Chunk [    ]) = liftI $ scan k rec_pos_ ref_pos
    scan k !rec_pos_ !ref_pos (Chunk (r:rs)) =
        case r of SNP {} -> let (_,!base) : (!qual,_) : _ = sort $ call (glf_lk r)
                            in ( if qual >= min_qual
                                 then lift $ enumPure1Chunk (ref ref_pos (rec_pos - ref_pos)) k
                                             >>= enumPure1Chunk [(base,qual)]
                                 else lift $ enumPure1Chunk (ref ref_pos (1 + rec_pos - ref_pos)) k )
                               >>= \k' -> scan k' rec_pos (1+rec_pos) (Chunk rs)

                  Indel {} | ins && iqual >= min_qual     -> lift (enumPure1Chunk (ref ref_pos (rec_pos + 1 - ref_pos)) k >>=
                                                                   enumPure1Chunk [ (b,iqual) | b <- S.unpack sq ]) >>= \k'' ->
                                                             scan k'' rec_pos ref_pos (Chunk rs)
                           | not ins && iqual >= min_qual -> lift (enumPure1Chunk (ref ref_pos (rec_pos - ref_pos)) k) >>= \k' ->
                                                             scan k' rec_pos (ref_pos + S.length sq) (Chunk (drop (S.length sq) rs))
                           | otherwise                    -> lift (enumPure1Chunk (ref ref_pos (rec_pos - ref_pos)) k) >>= \k' ->
                                                             scan k' rec_pos ref_pos (Chunk rs)
      where
        !rec_pos = rec_pos_ + glf_offset r
        (ins,sq) = if glf_lk_hom1 r > glf_lk_hom2 r
                   then (glf_is_ins2 r, glf_seq2 r) else (glf_is_ins1 r, glf_seq1 r)
        iqual = abs $ glf_lk_hom1 r - glf_lk_hom2 r


diploid_call, haploid_call :: [Int] -> [(Int, Char)]
diploid_call lks = zip lks "AMRWCSYGKT"
haploid_call lks = zip (map (lks !!) [0,4,7,9]) "ACGT"


type Formatter = String -> Enumeratee QSeq String IO ()

print_fasta :: Formatter
print_fasta name = eneeCheckIfDone (\k -> I.mapStream fst ><> toLines 60 $ k $ Chunk ('>' : name ++ "\n"))

print_fastq :: Formatter
print_fastq name = eneeCheckIfDone p'header
  where
    p'header k  = p'seq . k $ Chunk ('@' : name ++ "\n")
    p'seq it    = I.zip ((I.mapStream fst ><> toLines 60) it) (liftI $ coll [])
                  >>= \(it', qs) -> eneeCheckIfDone (p'sep qs) it'
    p'sep qs k  = lift $ (enumList (map S.unpack qs) >=> run) (toLines 60 . k $ Chunk "+\n")

    mkqual = chr . max 33 . min 126 . (+) 33 . fromIntegral
    coll !acc (EOF x) = lift (print $ length acc) >> idone (reverse acc) (EOF x)
    coll !acc (Chunk []) = liftI $ coll acc
    coll !acc (Chunk  c) = liftI . coll $! norm (S.pack (map (mkqual . snd) c)) acc

    -- ensure that we don't build many small ByteStrings
    norm !x [] = [x]
    norm !x (y:ys) | S.length x > S.length y = norm (y `S.append` x) ys
                   | otherwise               = x:y:ys


toLines :: Monad m => Int -> Enumeratee String String m r
toLines n = eneeCheckIfDone (\k -> I.isFinished >>= go k)
  where
    go k  True = return $ liftI k
    go k False = do s <- I.take n I.stream2list >>= lift . run
                    eneeCheckIfDone (\k1 -> toLines n . k1 $ Chunk "\n") . k $ Chunk s


readFasta :: L.ByteString -> [(L.ByteString, L.ByteString)]
readFasta = rd . dropWhile (not . isHeader) . L.lines
  where
    isHeader l = not (L.null l) && L.head l == '>'
    rd [] = []
    rd (l:ls) = let name = L.takeWhile (not . isSpace) $ L.drop 1 l
                    (sqs,rest) = break isHeader ls
                in (name, L.filter (`elem` "ACGTBDHVSWMKYRNU") $ L.concat sqs) : rd rest

subst_name :: String -> S.ByteString -> String
subst_name [] _ = []
subst_name ('%':'s':t) s = S.unpack s ++ subst_name t s
subst_name ('%':'%':t) s = '%' : subst_name t s
subst_name (t:ts) s = t : subst_name ts s

