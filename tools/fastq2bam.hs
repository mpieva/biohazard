{-# LANGUAGE BangPatterns, OverloadedStrings #-}
import Bio.Bam
import Bio.Bam.Evan ( removeWarts )
import Bio.Iteratee.ZLib
import Bio.Prelude
import System.Console.GetOpt
import System.IO

import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as S
import qualified Data.Vector.Generic as V

-- TODO:
-- - optional(!) GZip

data Opts = Opts { output :: BamMeta -> Iteratee [BamRec] IO ()
                 , inputs :: [Input]
                 , verbose :: Bool }

defaultOpts :: Opts
defaultOpts = Opts { output = protectTerm . pipeBamOutput
                   , inputs = []
                   , verbose = False }

data Input = Input { _read1 :: FilePath
                   ,  read2 :: Maybe FilePath
                   , index1 :: Maybe FilePath
                   , index2 :: Maybe FilePath }
  deriving Show

getopts :: [String] -> ([Opts -> IO Opts], [String], [String])
getopts = getOpt (ReturnInOrder add_read1) options
  where
    options =
        [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE"
        , Option "1" ["read-one"] (ReqArg add_read1 "FILE") "Parse FILE for anything"
        , Option "2" ["read-two"] (ReqArg add_read2 "FILE") "Parse FILE for second mates"
        , Option "I" ["index-one"] (ReqArg add_idx1 "FILE") "Parse FILE for first index"
        , Option "J" ["index-two"] (ReqArg add_idx2 "FILE") "Parse FILE for second index"
        , Option "v" ["verbose"] (NoArg set_verbose) "Print progress information"
        , Option "h?" ["help","usage"] (NoArg usage) "Print this helpful message" ]

    set_output "-" c = return $ c { output = pipeBamOutput }
    set_output  fn c = return $ c { output = writeBamFile fn }
    set_verbose    c = return $ c { verbose = True }

    add_read1 fn c = return $ c { inputs = Input fn Nothing Nothing Nothing : inputs c }
    add_read2 fn c = return $ c { inputs = at_head (\i -> i { read2  = Just fn }) (inputs c) }
    add_idx1  fn c = return $ c { inputs = at_head (\i -> i { index1 = Just fn }) (inputs c) }
    add_idx2  fn c = return $ c { inputs = at_head (\i -> i { index2 = Just fn }) (inputs c) }

    at_head f [    ] = [ f $ Input "-" Nothing Nothing Nothing ]
    at_head f (i:is) = f i : is

    usage _ = do pn <- getProgName
                 let t = "Usage: " ++ pn ++ " [OPTION...]\n" ++
                         "Reads multiple FastA or FastQ files and converts them to BAM.  See manpage for details."
                 hPutStrLn stderr $ usageInfo t options
                 exitSuccess


main :: IO ()
main = do (opts, [], errors) <- getopts `fmap` getArgs
          unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
          conf <- foldl (>>=) (return defaultOpts) opts
          pgm <- addPG Nothing

          let eff_inputs = if null (inputs conf) then [ Input "-" Nothing Nothing Nothing ] else inputs conf
          hPrint stderr $ eff_inputs

          foldr ((>=>) . enum_input) run (reverse eff_inputs) $
                joinI $ progress (verbose conf) $
                joinI $ mapChunks concatDuals $
                ilift liftIO $ output conf (pgm mempty)


type UpToTwo a = (a, Maybe a)

one :: a -> UpToTwo a
one a = (a, Nothing)

two :: a -> a -> UpToTwo a
two a b = (a, Just b)

mapU2 :: (a -> b) -> UpToTwo a -> UpToTwo b
mapU2 f (a,b) = (f a, fmap f b)

concatDuals :: [UpToTwo a] -> [a]
concatDuals ((a,Just  b):ds) = a : b : concatDuals ds
concatDuals ((a,Nothing):ds) = a : concatDuals ds
concatDuals [              ] = []

-- Enumerates a file.  Sequence and quality end up in b_seq and b_qual.
fromFastq :: (MonadIO m, MonadMask m) => FilePath -> Enumerator [BamRec] m a
fromFastq fp = enumAny fp $= enumInflateAny $= parseFastqCassava $= mapStream removeWarts
  where
    enumAny "-" = enumHandle defaultBufSize stdin
    enumAny  f  = enumFile defaultBufSize f

enum_input :: (MonadIO m, MonadMask m) => Input -> Enumerator [UpToTwo BamRec] m a
enum_input inp@(Input r1 mr2 mi1 mi2) o = do
    liftIO $ hPrint stderr inp
    (withIndex mi1 "XI" "YI" $ withIndex mi2 "XJ" "YJ" $
        case mr2 of Nothing -> fromFastq r1 $= mapStream one ; Just r2 -> enumDual r1 r2) o

-- Given an enumerator and maybe a filename, read index sequences from
-- the file and merge them with the numerator.
withIndex :: (MonadIO m, MonadMask m)
          => Maybe FilePath -> BamKey -> BamKey
          -> Enumerator [UpToTwo BamRec] m a -> Enumerator [UpToTwo BamRec] m a
withIndex Nothing      _    _ enum = enum
withIndex (Just fp) tagi tagq enum = mergeEnums enum (fromFastq fp) (convStream combine)
  where
    combine = do seqrecs <- lift headStream
                 idxrec  <- headStream
                 when (b_qname (fst seqrecs) /= b_qname idxrec) . error $
                        "read names do not match: " ++ shows (b_qname (fst seqrecs)) " & " ++ show (b_qname idxrec)

                 let idxseq  = S.pack $ map showNucleotides $ V.toList $ b_seq idxrec
                     idxqual = B.pack $ map   ((+33) . unQ) $ V.toList $ b_qual idxrec
                 return [ flip mapU2 seqrecs $
                        \r -> r { b_exts = (if B.null idxqual then id else insertE tagq (Text idxqual))
                                         $ insertE tagi (Text idxseq) $ b_exts r } ]

-- Enumerate dual files.  We read two FastQ files and match them up.  We
-- must make sure the names match, and we will flag everything as
-- 1st/2nd mate, no matter if the syntactic warts were present in the
-- files themselves.
enumDual :: (MonadIO m, MonadMask m)
         => FilePath -> FilePath -> Enumerator [UpToTwo BamRec] m a
enumDual f1 f2 = mergeEnums (fromFastq f1 $= mapStream one) (fromFastq f2) (convStream combine)
  where
    combine = do (firstMate, Nothing) <- lift headStream
                 secondMate           <- headStream

                 when (b_qname firstMate /= b_qname secondMate) . error $
                        "read names do not match: " ++ shows (b_qname firstMate) " & " ++ show (b_qname secondMate)

                 let qc = (b_flag firstMate .|. b_flag secondMate) .&. flagFailsQC
                     addx k = maybe id (updateE k) $ maybe (lookup k (b_exts secondMate)) Just $ lookup k (b_exts firstMate)
                     add_indexes = addx "XI" . addx "XJ" . addx "YI" . addx "YJ"

                 return [ two (firstMate  { b_flag = qc .|.  flagFirstMate .|. flagPaired .|. b_flag firstMate .&. complement flagSecondMate
                                          , b_exts = add_indexes $ b_exts firstMate })
                              (secondMate { b_flag = qc .|. flagSecondMate .|. flagPaired .|. b_flag secondMate .&. complement flagFirstMate
                                          , b_exts = add_indexes $ b_exts secondMate }) ]


progress :: MonadIO m => Bool -> Enumeratee [UpToTwo BamRec] [UpToTwo BamRec] m b
progress False = mapChunks id
progress True  = eneeCheckIfDonePass (icont . go 0 0)
  where
    go !_ !_ k (EOF         mx) = idone (liftI k) (EOF mx)
    go !l !n k (Chunk    [   ]) = liftI $ go l n k
    go !l !n k (Chunk as@(a:_)) = do
        let !n' = n + length as
            !nm = b_qname (fst a)
            !l' = l `max` S.length nm
        when (n `div` 2048 /= n' `div` 2048) $ liftIO $ do
            hPutStr stderr $ "\27[K" ++
                replicate (l' - S.length nm) ' '
                ++ S.unpack nm ++ ", "
                ++ shows n' " records processed\n"
            hFlush stderr
        eneeCheckIfDonePass (icont . go l' n') . k $ Chunk as


