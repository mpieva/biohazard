{-# LANGUAGE OverloadedStrings, BangPatterns, RecordWildCards #-}
-- Simple Mitochondrial Contamination Check on BAM files.
--
-- This is based on Ye Olde Contamination Check for the Neanderthal
-- genome; the method is the same (and will continue to not work on
-- modern humans), but simplified and sanitized.  Differences from
-- before:
--
-- * We use the alignment from the BAM file as is.  Earlier we would
--   have created *two* new alignments.  That is silly, however.  Two
--   new alignments should not be followed by bean counting, but by an
--   attempt to genotype both the sample and the contaminant.
--
-- * Before, the sample and contaminant sequences were fixed.  Now we
--   instead input a list of the diagnostic positions.  Instead of an
--   explicit list, the two sequences can still be used, or only the
--   contaminant can be supplied while the sample is genotype-called.


-- TODO
--
-- (1) Given a list of diagnostic positions, implement the contamination
--     check.  Structure of the code can be stolen from ccheck in the
--     mia package.
--
--     - What do we do about wrapped alignments?  Mia has f/b/a labels,
--       BAM doesn't.  We can see if it overhangs, though.
--
-- (2) Given a high-coverage sample, genotype call it and derive the
--     diagnostic positions.
--
--     - This method needs some definition of the contaminant consensus
--       thingy.
--
-- (3) Given a `correct' sample sequence, align it to the reference and
--     derive diagnostic positions from that.
--
--     - Needs the same description of the contaminant thingy.
--
-- (4) Consider Read Groups.
--
--     - One result per read group (or maybe per library, alternatively
--       per file) should be produced.
--     - The "aDNA" setting should be determined from either the @RG
--       header or from an external source.


import Bio.Base
import Bio.Bam
import Bio.Iteratee
import Control.Applicative
import Control.Monad
import Data.Bits
import Data.Monoid
import Data.List
import Numeric
import System.Console.GetOpt
import System.Environment
import System.Exit
import System.IO

import qualified Data.HashMap.Strict        as HM
import qualified Data.IntMap                as IM
import qualified Data.ByteString.Char8      as S

data Conf = Conf {
        conf_adna :: Adna,
        conf_verbosity :: Int,
        conf_header :: HeaderFn,
        conf_output :: OutputFn,
        conf_shoot_foot :: Bool,
        conf_dp_list :: DpList
    }

options :: [OptDescr (Conf -> IO Conf)]
options = public_options ++ hidden_options
  where
    public_options = [
        Option "a" ["ancient","dsprot"] (NoArg (set_adna ancientDNAds)) "Treat DNA as ancient, double strand protocol",
        Option "s" ["ssprot"]           (NoArg (set_adna ancientDNAss)) "Treat DNA as ancient, single strand protocol",
        Option ""  ["fresh"]            (NoArg (set_adna     freshDNA)) "Treat DNA as fresh (not ancient)",
        Option "T" ["table"]            (NoArg        set_output_table) "Print output in table form",

        Option "v" ["verbose"]          (NoArg    (mod_verbosity succ)) "Produce more debug output",
        Option "q" ["quiet"]            (NoArg    (mod_verbosity pred)) "Produce less debug output",
        Option "h?" ["help","usage"]    (NoArg                   usage) "Print this message and exit"
      ]

    hidden_options = [
        Option ""  ["shoot","foot"]     (NoArg set_shoot_foot) []
      ]

    usage _ = do pn <- getProgName
                 hPutStrLn stderr $ usageInfo ("Usage: " ++ pn ++ " [OPTION...] [Bam-File...]") public_options
                 exitSuccess

    set_shoot_foot   c = return $ c { conf_shoot_foot = True }
    set_adna       a c = return $ c { conf_adna = a }
    set_output_table c = return $ c { conf_output = show_result_table, conf_header = header_table }
    mod_verbosity  f c = return $ c { conf_verbosity = f (conf_verbosity c) }


conf0 :: IO Conf
conf0 = return $ Conf { conf_adna = freshDNA
                      , conf_verbosity = 1
                      , conf_header = ""
                      , conf_output = show_result_plain
                      , conf_shoot_foot = False
                      , conf_dp_list = error "no diagnostic positions defined"
                      }

{- Old options... may or may not be of much use.

struct option longopts[] = {
	{ "reference", required_argument, 0, 'r' },
	{ "transversions", no_argument, 0, 't' },
	{ "span", required_argument, 0, 's' },
	{ "maxd", required_argument, 0, 'd' },
} ;

void usage( const char* pname )
{
		"Reads a maln file and tries to quantify contained contamination.\n"
		"Options:\n"
		"  -r, --reference FILE     FASTA file with the likely contaminant (default: builtin mt311)\n"
		"  -t, --transversions      Treat only transversions as diagnostic\n"
		"  -s, --span M-N           Look only at range from M to N\n"
		"  -n, --numpos N           Require N diagnostic sites in a single read (default: 1)\n"
}
-}

-- | A list of diagnostic positions.  We drop the failed idea of
-- "weakly diagnostic positions".  We also work in the coordinate system
-- of the reference.  Therefore, a diagnostic position is defined by
-- position, allele in the clean sample and allele in the contaminant.

data Dp = Dp { dp_clean_allele :: !Nucleotide
             , dp_dirty_allere :: !Nucleotide }
  deriving Show

type DpList = IM.IntMap Dp

show_dp_list :: DpList -> ShowS
show_dp_list = flip IM.foldrWithKey id $ \pos (Dp cln drt) k ->
    (:) '<' . shows pos . (:) ':' . shows drt .
    (:) ',' . shows cln . (++) ">, " . k


-- | Reads are classified into one of these.
data Klass = Unknown | Clean | Dirty | Conflict | Nonsense
  deriving (Ord, Eq, Enum, Bounded, Show)

instance Monoid Klass where
    mempty = Unknown
    Clean `mappend` Dirty = Conflict
    Dirty `mappend` Clean = Conflict
    x `mappend` y = if x < y then y else x

newtype Summary = Summary (IM.IntMap Int)

sum_count :: Klass -> Summary -> Summary
sum_count kl (Summary m) = Summary $ IM.insertWith (+) (fromEnum kl) 1 m

sum_get :: Klass -> Summary -> Int
sum_get kl (Summary m) = IM.findWithDefault 0 (fromEnum kl) m


-- | Determines what an allele could come from.  Does not take
-- port-mortem modifications into account.
classify :: Dp -> Nucleotide -> Klass
classify (Dp cln drt) nuc
    | maybe_clean && maybe_dirty = Unknown
    | maybe_clean                = Clean
    |                maybe_dirty = Dirty
    | otherwise                  = Nonsense
  where
    maybe_clean = unN cln .&. unN nuc /= 0
    maybe_dirty = unN drt .&. unN nuc /= 0


-- | We deal with aDNA by transforming a base into all the bases it
-- could have been.  So the configuration is simply the transformation
-- function.
type Adna = Nucleotide -> Nucleotide

-- | Fresh DNA: no transformation.
freshDNA :: Adna
freshDNA = id

-- | Ancient DNA, single strand protocol.  Deamination can turn C into T
-- only.
ancientDNAss :: Adna
ancientDNAss = N . app . unN
  where app x = if x .&. unN nucT /= 0 then x .|. unN nucC else x

-- | Ancient DNA, double strand protocol.  Deamination can turn C into T
-- and G into A.
ancientDNAds :: Adna
ancientDNAds = N . app1 . app2 . unN
  where app1 x = if x .&. unN nucT /= 0 then x .|. unN nucC else x
        app2 x = if x .&. unN nucA /= 0 then x .|. unN nucG else x


-- | Classifying a read.  In an ideal world, we'd be looking at a single
-- read mapped in one piece.  Instead, we may be looking at half a mate
-- pair or even a single read mapped inconveniently across the origin.
--
-- We will be reading a BAM stream.  All reads with the same name (there
-- maybe 1..4, assuming no major breakage) need to be processed
-- together.  We'll isolate that here:  our input stream consists of
-- reads that all have the same qname.  Results in exactly one 'Klass'.
-- We will ignore mate pairs that are improperly mapped or filtered.
--
-- May need more options.  Note that application of the aDNA function
-- depends on the strandedness of the alignment.  FIXME
--
-- This is the only place where counting of votes was used before, and
-- only for debugging purposes.  Everything that was either dirty or
-- clean (but not both) counted as a vote.

classify_read_set :: Monad m => DpList -> Adna -> Iteratee [BamRaw] m Klass
classify_read_set = undefined

-- | Classifying a stream.  We create a map from read name to iteratee.
-- New names are inserted, known names fed to stored iteratees.
-- ``Done'' iteratees are disposed of immediately.

classify_stream :: Monad m => DpList -> Adna -> Iteratee [BamRaw] m Summary
classify_stream dps adna = foldStreamM classify_read (Summary IM.empty, HM.empty) >>= lift . finish
  where
    classify0 = classify_read_set dps adna

    classify_read (summary, iters) rd = do
        let it = HM.lookupDefault classify0 (br_qname rd) iters
        (isdone, it') <- enumPure1Chunk [rd] it >>= enumCheckIfDone
        if isdone then do cl <- run it'
                          return (sum_count cl summary, HM.delete (br_qname rd) iters)
                  else return (summary, HM.insert (br_qname rd) it' iters)

    finish (summary, iters) = foldM (\s it -> flip sum_count s `liftM` run it) summary $ HM.elems iters


{- Missing from the output right now:

 * filename (library would be better)
 * alignment distance (only useful if DPs are derived from alignment)
 * number of difference (likewise)
 * number of DPs
 * number of DPs which are transversions
-}

result_labels :: [ String ]
result_labels = [ "unclassified", "clean", "polluting", "conflicting", "nonsensical", "LB", "ML", "UB" ]

type HeaderFn = String
type OutputFn = Summary -> Maybe [Double] -> String

show_result_plain :: OutputFn
show_result_plain summary ests = unlines $ zipWith fmt result_labels [minBound..maxBound] ++ [[]]
  where
    labellen = (+) 2 . maximum . map length $ zipWith const result_labels [minBound..maxBound::Klass]
    pad n s  = replicate (n - length s) ' ' ++ s

    fmt lbl kl = pad labellen lbl ++ " fragments: " ++ show (sum_get kl summary) ++
                 if kl == Dirty then maybe [] fmt_ests ests else []

    fmt_ests [lb,ml,ub] = " (" ++ showFFloat (Just 1) lb " .. "
                               ++ showFFloat (Just 1) ml " .. "
                               ++ showFFloat (Just 1) ub "%)"

header_table :: HeaderFn
header_table = intercalate "\t" result_labels

show_result_table :: OutputFn
show_result_table summary ests = intercalate "\t" $
    [ show $ sum_get kl summary | kl <- [minBound..maxBound] ] ++
    maybe (replicate 3 "N/A") (map (\x -> showFFloat (Just 1) x [])) ests


show_result_with :: (Summary -> Maybe [Double] -> a) -> Summary -> a
show_result_with f summary = f summary (if nn /= 0 then Just [lb,ml,ub] else Nothing)
  where
    z = 1.96   -- this is Z_{0.975}, giving a 95% confidence interval
    k =     fromIntegral (sum_get Dirty summary)
    n = k + fromIntegral (sum_get Clean summary)
    nn = sum_get Dirty summary + sum_get Clean summary

    p_ = k / n
    c = p_ + 0.5 * z * z / n
    w = z * sqrt( p_ * (1-p_) / n + 0.25 * z * z / (n*n) )
    d = 1 + z * z / n

    lb = max  0  $ 100 * (c-w) / d    -- lower bound of CI
    ml =           100 * p_           -- ML estimate
    ub = min 100 $ 100 * (c+w) / d    -- upper bound of CI


-- The following is old 'ccheck'... for reference and guidance.


{-
/*
 * Contamination Checker.  Outline:
 *
 * - read the human reference (concsensus of contaminants); this will
 *   contain ambiguity codes
 * - read maln file, including assembly and assembled reads
 * - align contaminant-consensus and assembly globally
 *   This uses Myers' O(nd) aligner, for it grasps ambiguity codes and
 *   runs fast enough, in little memory, for long, but similar
 *   sequences.
 * - find "strongly diagnostic positions", positions where ass and con
 *   are incompatible, and "weakly diagnostic positions", positions
 *   where ass and con are not always equal
 * - for every "end" fragment: store it  and later join with its other
 *   half to give an effectively "full" fragment
 * - for every "full" fragment: if it crosses at least one (strongly or
 *   weakly) diagnostic position, cut out that range from ref and align
 *   to it globally using the mia aligner
 * - pass 1: for every weakly diagnostic position where the bases agree,
 *   store whether a contaminant was discovered, and if so, turn them
 *   into "actually diagnostic positions".
 * - pass 2: for every (strongly or actually) diagnostic position where
 *   the bases agree, classify it, then classify the fragment
 *   (conflicting, uninformative, contaminant, endogenous)
 * - produce a summary
 *
 * Notable features:
 * - operates sensibly on aDNA
 * - has sensible commandline and doesn't make too much noise in operation
 * - optionally considers only certain diagnostic positions
 *   (tranversions only and/or some region only)
 * - new consensus sequence has other letters besides N
 */

// Everything that differs is weakly diagnostic, unless it's a gap.
// Note that this mean that Ns are usually weakly diagnostic.
bool is_diagnostic( char aln1, char aln2 )
{
	return aln1 != '-' && aln2 != '-' && toupper(aln1) != toupper(aln2) ;
}

// Interesting question... given ambiguity codes, what's a transversion?
// One way to put it:  anything that is incompatible with all four
// transitions.  Needs a different implementation.
bool is_transversion( char a, char b )
{
	char u = a & ~32 ;
	char v = b & ~32 ;
	switch( u )
	{
		case 'A': return v != 'G' ;
		case 'C': return v != 'T' ;
		case 'G': return v != 'A' ;
		case 'T':
		case 'U': return v != 'C' ;
		default: return false ;
	}
}


dp_list mk_dp_list( const char* aln1, const char* aln2, int span_from, int span_to )
{
	dp_list l ;
    int index = 0 ;
    while( index != span_from && *aln1 && *aln2 )
    {
		if( *aln2 != '-' ) ++index ;
		++aln1 ;
		++aln2 ;
    }
	while( index != span_to && *aln1 && *aln2 )
	{
		if( is_diagnostic( *aln1, *aln2 ) ) {
            l[index].consensus = *aln1 ;
            l[index].assembly = *aln2 ;
        }
		if( *aln2 != '-' ) ++index ;
		++aln1 ;
		++aln2 ;
	}
	return l ;
}
-}

-- We won't keep this.  Mt311 should be stored as half a Dp list.
-- extern       char mt311_sequence[] ;
-- extern const int  mt311_sequence_size ;

main :: IO ()
main = do
    (opts, files, errors) <- getOpt Permute options <$> getArgs
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Conf{..} <- foldl (>>=) conf0 opts

{-
	bool transversions = false ;
	int min_diag_posns = 1 ;
	int maxd = 0 ;
	int span_from = 0, span_to = INT_MAX ;

	int opt ;
	do {
		opt = getopt_long( argc, argv, "r:avhts:d:n:MfTF", longopts, 0 ) ;
		switch( opt )
		{
			case 'r': read_fasta_ref( &hum_ref, optarg ) ; break ;
			case 't': transversions = true ; break ;
			case 's': sscanf( optarg, "%u-%u", &span_from, &span_to ) ; if( span_from ) span_from-- ; break ;
			case 'n': min_diag_posns = atoi( optarg ) ; break ;
			case 'd': maxd = atoi( optarg ) ; break ;
		}
	} while( opt != -1 ) ;
-}

    when (IM.size conf_dp_list < 40 && not conf_shoot_foot) $ do
        hPutStrLn stderr $
            "\n *** Low number (" ++ shows (IM.size conf_dp_list) ") of diagnostic positions found.\n\
              \ *** I will stop now for your own safety.\n\
              \ *** If you are sure you want to shoot yourself\n\
              \ *** in the foot, read the man page to learn\n\
              \ *** how to lift this restriction.\n\n"
        exitFailure

    -- TODO  We will usually want to seek to the mitochondrion, which
    -- doesn't work with the simple 'mergeInputs' invocation.
    r <- mergeInputs combineCoordinates files >=> run $ \hdr ->
            classify_stream conf_dp_list conf_adna

    putStrLn $ unlines $ conf_header : show_result_with conf_output r : []

        {-
        if( mktable ) {
            fputs( infile.c_str(), stdout ) ;
            putchar( '\t' ) ;
        }
        else {
            puts( infile.c_str() ) ;
            putchar( '\n' ) ;
        }
        -}

        -- if( !maxd ) maxd = max( strlen(hum_ref.seq), strlen(maln->ref->seq) ) / 10 ;
--         char *aln_con = (char*)malloc( strlen(hum_ref.seq) + maxd + 2 ) ;
  --       char *aln_ass = (char*)malloc( strlen(maln->ref->seq) + maxd + 2 ) ;
    --     unsigned d = myers_diff( hum_ref.seq, myers_align_globally, maln->ref->seq, maxd, aln_con, aln_ass ) ;

        {-
        if( d == UINT_MAX ) {
            fprintf( stderr, "\n *** Could not align references with up to %d mismatches.\n"
                             " *** This is usually a sign of trouble, but\n"
                             " *** IF AND ONLY IF YOU KNOW WHAT YOU ARE DOING, you can\n"
                             " *** try the -d N option with N > %d.\n\n", maxd, maxd ) ;
            return 1 ;
        }
        if( mktable ) printf( "%d\t", d ) ;
        else printf( "  %d alignment distance between reference and assembly.\n", d ) ;

        if( verbose >= 6 ) print_aln( aln_con, aln_ass ) ;

        dp_list l = mk_dp_list( aln_con, aln_ass, span_from, span_to ) ;
        if( mktable ) printf( "%u\t", (unsigned)l.size() ) ;
        else printf( "  %u total differences between reference and assembly.\n", (unsigned)l.size() ) ;

        int num_strong = 0 ;
        for( dp_list::const_iterator i = l.begin() ; i != l.end() ; ++i )
            if( i->second.strength > weak ) ++num_strong ;
        if( mktable ) printf( "%d\t", (int)l.size() ) ;
        else {
            printf( "  %d diagnostic positions", (int)l.size() ) ;
            if( span_from != 0 || span_to != INT_MAX )
                printf( " in range [%d,%d)", span_from, span_to ) ;
            printf( ", %d of which are strongly diagnostic.\n", num_strong ) ;
        }

        if( verbose >= 3 ) {
            print_dp_list( stderr, l.begin(), l.end(), '\n', 0 ) ;
            print_dp_list( stderr, l.begin(), l.end(), '\n', 1 ) ;
        }

-}

        {-
        if( verbose >= 2 ) fputs( "Pass one: finding actually diagnostic positions.\n", stderr ) ;
        for( const AlnSeqP *s = maln->AlnSeqArray ; s != maln->AlnSeqArray + maln->num_aln_seqs ; ++s )
        {
            fixup_name( *s ) ;

            std::string the_ass( maln->ref->seq + (*s)->start, (*s)->end - (*s)->start + 1 ) ;
            // are we overlapping anything at all?
            std::pair< dp_list::const_iterator, dp_list::const_iterator > p =
                overlapped_diagnostic_positions( l, *s ) ;

            if( verbose >= 3 )
            {
                fprintf( stderr, "%s/%c:\n  %d potentially diagnostic positions",
                         (*s)->id, (*s)->segment, (int)std::distance( p.first, p.second ) ) ;
                if( verbose >= 4 )
                {
                    putc( ':', stderr ) ; putc( ' ', stderr ) ;
                    print_dp_list( stderr, p.first, p.second, 0 ) ;
                }
                fprintf( stderr, "; range:  %d..%d\n", (*s)->start, (*s)->end ) ;
            }
-}


        {-
            int t = 0 ;
            for( dp_list::const_iterator i = l.begin() ; i != l.end() ; ++i )
                if( is_transversion( i->second.consensus, i->second.assembly ) ) ++t ;
            if( mktable ) printf( "%d\t%d\t", t, num_strong ) ;
            else {
                printf( "  %d effectively diagnostic positions", (int)l.size() ) ;
                if( span_from != 0 || span_to != INT_MAX )
                    printf( " in range [%d,%d)", span_from, span_to ) ;
                printf( ", %d of which are transversions.\n\n", t ) ;
            }
        if( verbose >= 3 ) print_dp_list( stderr, l.begin(), l.end(), '\n' ) ;

        std::deque< cached_pwaln >::const_iterator cpwaln = cached_pwalns.begin() ;
        for( const AlnSeqP *s = maln->AlnSeqArray ; s != maln->AlnSeqArray + maln->num_aln_seqs ; ++s, ++cpwaln )
        {
            whatsit klass = unknown ;
            int votes = 0, votes2 = 0 ;

            std::string the_ass( maln->ref->seq + (*s)->start, (*s)->end - (*s)->start + 1 ) ;
            // enough overlap?  (we only have _actually_ diagnostic positions now)
            std::pair< dp_list::const_iterator, dp_list::const_iterator > p =
                overlapped_diagnostic_positions( l, *s ) ;
            if( std::distance( p.first, p.second ) < min_diag_posns )
            {
                if( verbose >= 3 ) {
                    fputs( (*s)->id, stderr ) ;
                    putc( '/', stderr ) ;
                    putc( (*s)->segment, stderr ) ;
                    fputs( ": no diagnostic positions\n", stderr ) ;
                }
            }
            else
            {
                if( verbose >= 3 )
                {
                    fprintf( stderr, "%s/%c: %d diagnostic positions", (*s)->id, (*s)->segment, (int)std::distance( p.first, p.second ) ) ;
                    if( verbose >= 4 )
                    {
                        putc( ':', stderr ) ; putc( ' ', stderr ) ;
                        print_dp_list( stderr, p.first, p.second, 0 ) ;
                    }
                    fprintf( stderr, "; range:  %d..%d\n", (*s)->start, (*s)->end ) ;
                }

                // Hmm, all this iterator business is somewhat lacking...
                char *paln1 = aln_con, *paln2 = aln_ass ;
                int ass_pos = 0 ;
                while( ass_pos != (*s)->start && *paln1 && *paln2 )
                {
                    if( *paln2 != '-' ) ass_pos++ ;
                    ++paln1 ;
                    ++paln2 ;
                }

                char *in_ass = maln->ref->seq + (*s)->start ;
                char *in_frag_v_ass = (*s)->seq ;
                std::string::const_iterator in_frag_v_ref = cpwaln->frag_seq.begin() ;

                std::string lifted = lift_over( aln_con, aln_ass, (*s)->start, (*s)->end + 1 ) ;
                std::string in_ref = lifted.substr( 0, cpwaln->start ) ;
                in_ref.append( cpwaln->ref_seq ) ;

                while( ass_pos != (*s)->end +1 && *paln1 && *paln2 && !in_ref.empty() && *in_ass && *in_frag_v_ass && *in_frag_v_ref )
                {
                    if( *paln1 != '-' ) {
                        do {
                            in_ref=in_ref.substr(1) ;
                            in_frag_v_ref++ ;
                        } while( in_ref[0] == '-' ) ;
                    }
                    if( *paln2 != '-' ) {
                        ass_pos++ ;
                        do {
                            in_ass++ ;
                            in_frag_v_ass++ ;
                        } while( *in_ass == '-' ) ;
                    }
                    ++paln1 ;
                    ++paln2 ;
                }
                if( verbose >= 4 ) putc( '\n', stderr ) ;
            }
        }
    }
}
        -}

