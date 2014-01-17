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

import Bio.Base
import Bio.Bam
import Bio.Iteratee

import qualified Data.HashMap.Strict as HM
import qualified Data.IntMap.Strict as IM
import qualified Data.ByteString as B

data Conf = Conf {}

options :: [OptDescr (Conf -> IO Conf)]
options = []

{- Old options... may or may not be of much use.


struct option longopts[] = {
	{ "reference", required_argument, 0, 'r' },
	{ "ancient", no_argument, 0, 'a' },
	{ "verbose", no_argument, 0, 'v' },
	{ "help",    no_argument, 0, 'h' },
	{ "transversions", no_argument, 0, 't' },
	{ "span", required_argument, 0, 's' },
	{ "maxd", required_argument, 0, 'd' },
    { "table", no_argument, 0, 'T' },
    { "shoot", no_argument, 0, 'F' },
    { "foot", no_argument, 0, 'F' },
	{ 0,0,0,0 }
} ;

void usage( const char* pname )
{
	fputs( "Usage: ", stdout ) ;
	fputs( pname, stdout ) ;
	fputs( " [-r <ref.fa>] [-a] [-t] [-s M-N] [-v] <aln.maln> \n\n"
		"Reads a maln file and tries to quantify contained contamination.\n"
		"Options:\n"
		"  -r, --reference FILE     FASTA file with the likely contaminant (default: builtin mt311)\n"
		"  -a, --ancient            Treat DNA as ancient (i.e. likely deaminated)\n"
		"  -t, --transversions      Treat only transversions as diagnostic\n"
		"  -s, --span M-N           Look only at range from M to N\n"
		"  -n, --numpos N           Require N diagnostic sites in a single read (default: 1)\n"
        "  -T, --table              Output as tables (easier for scripts, harder on the eyes)\n"
		"  -v, --verbose            Increase verbosity level (can be repeated)\n"
		"  -h, --help               Print this help message\n\n", stdout ) ;
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
show_dp_list = foldMapWithKey $ \pos (Dp cln drt) ->
    (:) '<' . shows pos . (:) ':' . shows drt .
    (:) ',' . shows cln . (++) ">, "


-- | Reads are classified into one of these.
data Klass = Unknown | Clean | Dirty | Conflict | Nonsense
  deriving (Ord, Eq, Enum, Show)

instance Monoid Klass where
    mempty = Unknown
    Clean `mappend` Dirty = Conflict
    Dirty `mappend` Clean = Conflict
    x `mappend` y = if x < y then y else x

newtype Summary = Summary (IM.IntMap Int)

sum_count :: Klass -> Summary -> Summary
sum_count kl (Summary m) = Summary $ IM.insertWith (+) (fromEnum kl) 1 m

-- | Determines what an allele could come from.  Does not take
-- port-mortem modifications into account.
classify :: Dp -> Nucleotide -> Klass
classify (Dp cln drt) nuc
    | maybe_cln && maybe_dirty = Unknown
    | maybe_cln                = Clean
    |              maybe_dirty = Dirty
    | otherwise                = Nonsense
  where
    maybe_clean = unN cln .&. unN nuc /= 0
    maybe_dirty = unN drt .&. unN nuc /= 0


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
-- Needs options.  At least the Dp list, probably an aDNA kind of
-- setting, too.
--
-- This is the only place where counting of votes was used before, and
-- only for debugging purposes.  Everything that was either dirty or
-- clean (but not both) counted as a vote.

classify_read_set :: Monad m => Iteratee [BamRaw] m Klass
classify_read_set = undefined

-- | Classifying a stream.  We create a map from read name to iteratee.
-- New names are inserted, known names fed to stored iteratees.
-- ``Done'' iteratees are disposed of immediately.

classify_stream :: Monad m => Iteratee [BamRaw] m Summary
classify_stream = foldStreamM classify_read (Summary IM.empty, HM.empty) >>= lift . finish
  where
    classify0 = classify_read_set

    classify_read (summary, iters) rd = do
        let it = HM.lookupDefault classify0 (br_qname rd) iters
        (isdone, it') <- lift $ enumPure1Chunk [rd] it >>= enumCheckIfDone
        if isdone then do cl <- lift $ run it'
                          return (sum_count summary, HM.delete (br_qname rd) iters)
                  else return (summary, HM.insert (br_qname rd) it' iters)

    finish (summary, iters) = foldM (\s it -> sum_count s <$> run it) summary $ HM.elems iters


result_labels :: [ S.ByteString ]
result_labels = [ "unclassified", "clean", "polluting", "conflicting", "nonsensical", "LB", "ML", "UB" ]

show_result_plain :: Summary -> Maybe [Double] -> S.ByteString
show_result_plain summary ests = S.unlines $ flip (:) S.empty $
    zipWith fmt labels [ sum_get kl summary | kl <- [minBound..maxBound] ] ++ [S.empty]
  where
    labellen = (+) 2 . maximum . map S.length . zipWith snd [minBound..maxBound::Klass] $ result_labels

    fmt lbl kl = pad labellen lbl ++ " fragments: " ++ show (sum_get kl summary) ++
                 if kl == Dirty then maybe [] fmt_ests ests else []

    fmt_ests [lb,ml,ub] = " (" ++ showFFloat (Just 1) lb " .. "
                               ++ showFFloat (Just 1) ml " .. "
                               ++ showFFloat (Just 1) ub "%)"

show_result_table :: Summary -> Maybe [Double] -> String
show_result_table summary ests = intercalate "\t" $
    [ show $ sum_get kl summary | kl <- (minBound, maxBound) ] ++
    maybe (replicate 3 "N/A") (map (\x -> showFFloat (Just 1) x []))


show_result_with :: (Summary -> Maybe [Double] -> a) -> Summary -> a
show_result_with k summary = k summary (if nn != 0 then Just [lb,ml,ub] else Nothing)
  where
    z = 1.96   -- this is Z_{0.975}, giving a 95% confidence interval
    k =     fromIntegral $ sum_get Dirty summary
    n = k + fromIntegral $ sum_get Clean summary
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

void print_aln( const char* aln1, const char* aln2 )
{
	int p ;
	const char* a ;
	while( *aln1 && *aln2 ) {
		for( p = 0, a = aln1 ; *a && p != 72 ; ++p ) putc( *a++, stderr ) ;
		putc( '\n', stderr ) ;

		for( p = 0, a = aln2 ; *a && p != 72 ; ++p ) putc( *a++, stderr ) ;
		putc( '\n', stderr ) ;

		for( p = 0 ; *aln1 && *aln2 && p != 72 ; ++p )
			putc( *aln1++ == *aln2++ ? '*' : ' ', stderr ) ;
		putc( '\n', stderr ) ;
		putc( '\n', stderr ) ;
	}
}


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

std::pair< dp_list::const_iterator, dp_list::const_iterator >
overlapped_diagnostic_positions( const dp_list& l, const AlnSeqP s )
{
	dp_list::const_iterator left  = l.lower_bound( s->start ) ;
	dp_list::const_iterator right = l.lower_bound( s->end + 1 ) ;
	return std::make_pair( left, right ) ;
}

// XXX: linear scan --> O(n)
// This could be faster (O(log n)) if precompiled into some sort of index.
std::string lift_over( const char* aln1, const char* aln2, int s, int e )
{
	std::string r ;
	int p ;
	for( p = 0 ; p < e && *aln1 && *aln2 ; ++aln1, ++aln2 )
	{
		if( *aln1 != '-' && p >= s ) r.push_back( *aln1 ) ;
		if( *aln2 != '-' ) ++p ;
	}
	return r ;
}

bool consistent( bool adna, char x, char y )
{
	char x_ = x == 'G' ? 'R' : x == 'C' ? 'Y' :
	          x == 'g' ? 'r' : x == 'c' ? 'y' : x ;
	return x == '-' || y == '-' || (char_to_bitmap( adna ? x_ : x ) & char_to_bitmap(y)) != 0 ;
}

// We won't keep this.  Mt311 should be stored a half a Dp list.
extern       char mt311_sequence[] ;
extern const int  mt311_sequence_size ;


int main( int argc, char * const argv[] )
{
	bool adna = false ;
	bool transversions = false ;
    bool be_clever = true ;
    bool mktable = false ;
    bool really = false ;
	int min_diag_posns = 1 ;
	int verbose = 0 ;
	int maxd = 0 ;
	int span_from = 0, span_to = INT_MAX ;

	if( argc == 0 ) { usage( argv[0] ) ; return 0 ; }

	int opt ;
	do {
		opt = getopt_long( argc, argv, "r:avhts:d:n:MfTF", longopts, 0 ) ;
		switch( opt )
		{
			case 'r':
				read_fasta_ref( &hum_ref, optarg ) ;
				break ;
			case 'a':
				adna = true ;
				break ;
			case 'v':
				++verbose ;
				break ;
			case ':':
				fputs( "missing option argument\n", stderr ) ;
				break ;
			case '?':
				fputs( "unknown option\n", stderr ) ;
				break ;
			case 'h':
				usage( argv[0] ) ;
				return 1 ;
			case 't':
				transversions = true ;
				break ;
			case 's':
				sscanf( optarg, "%u-%u", &span_from, &span_to ) ;
				if( span_from ) span_from-- ;
				break ;
			case 'n':
				min_diag_posns = atoi( optarg ) ;
				break ;
			case 'd':
				maxd = atoi( optarg ) ;
				break ;
            case 'M':
                break ;
            case 'f':
                be_clever = false ;
                break ;
            case 'T':
                mktable = true ;
                break ;
            case 'F':
                really = true;
                break ;
		}
	} while( opt != -1 ) ;

	if( optind == argc ) { usage( argv[0] ) ; return 1 ; }
	bool hum_ref_ok = sanity_check_sequence( hum_ref.seq ) ;
	if( !hum_ref_ok ) fputs( "FUBAR'ed FastA file: contaminant sequence contains gap symbols.\n", stderr ) ;

	if( !hum_ref.rcseq ) make_reverse_complement( &hum_ref ) ;

    if( mktable ) {
        fputs( "#Filename\tAln.dist\t#diff\t#weak\t#tv", stdout ) ;
        for( int i =0 ; i != 2 ; ++i ) {
            fputs( i ? "\t#eff" : "\t#strong", stdout ) ;
            for( int klass = 0 ; klass != sizeof(label)/sizeof(label[0]) ; ++klass )
            {
                putchar( '\t' ) ;
                fputs( label[klass], stdout ) ;
                if( i ) putchar( '\'' ) ;
            }
        }
        putchar( '\n' ) ;
    }

    for( ; optind != argc ; ++optind )
    {
        int summary[ maxwhatsits ] = {0} ;
        int summary2[ maxwhatsits ] = {0} ;

        std::string infile( be_clever ? find_maln( argv[optind] ) : argv[optind] ) ;
        if( mktable ) {
            fputs( infile.c_str(), stdout ) ;
            putchar( '\t' ) ;
        }
        else {
            puts( infile.c_str() ) ;
            putchar( '\n' ) ;
        }
        MapAlignmentP maln = read_ma( infile.c_str() ) ;
        PSSMP submat = maln->fpsm ;

        bool maln_ref_ok = sanity_check_sequence( maln->ref->seq ) ;
        if( !maln_ref_ok ) fputs( "FUBAR'ed maln file: consensus sequence contains gap symbols.\n", stderr ) ;
        if( !hum_ref_ok || !maln_ref_ok ) {
            fputs( "Problem might exist between keyboard and chair.  I give up.\n", stderr ) ;
            return 1 ;
        }

        if( !maxd ) maxd = max( strlen(hum_ref.seq), strlen(maln->ref->seq) ) / 10 ;
        char *aln_con = (char*)malloc( strlen(hum_ref.seq) + maxd + 2 ) ;
        char *aln_ass = (char*)malloc( strlen(maln->ref->seq) + maxd + 2 ) ;
        unsigned d = myers_diff( hum_ref.seq, myers_align_globally, maln->ref->seq, maxd, aln_con, aln_ass ) ;

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

        if( num_strong < 40 && !really ) {
            fprintf( stderr, "\n *** Low number (%d) of diagnostic positions found.\n"
                             " *** I will stop now for your own safety.\n"
                             " *** If you are sure you want to shoot yourself\n"
                             " *** in the foot, read the man page to learn\n"
                             " *** how to lift this restriction.\n\n", num_strong ) ;
            return 1 ;
        }

        typedef std::map< std::string, std::pair< whatsit, int > > Bfrags ;
        Bfrags bfrags, bfrags2 ;
        std::deque< cached_pwaln > cached_pwalns ;

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

            // reconstruct read and reference sequences, align them
            std::string the_read ;
            for( char *nt = (*s)->seq, **ins = (*s)->ins ; *nt ; ++nt, ++ins )
            {
                if( *nt != '-' ) the_read.push_back( *nt ) ;
                if( *ins ) the_read.append( *ins ) ;
            }
            std::string lifted = lift_over( aln_con, aln_ass, (*s)->start, (*s)->end + 2 ) ;

            if( verbose >= 5 )
            {
                fprintf( stderr, "\nraw read: %s\nlifted:   %s\nassembly: %s\n\n"
                        "aln.read: %s\naln.assm: %s\nmatches:  ",
                        the_read.c_str(), lifted.c_str(), the_ass.c_str(),
                        (*s)->seq, the_ass.c_str() ) ;
                std::string::const_iterator b = the_ass.begin(), e = the_ass.end() ;
                const char* pc = (*s)->seq ;
                while( b != e && *pc ) putc( *b++ == *pc++ ? '*' : ' ', stderr ) ;
            }

            int size = std::max( lifted.size(), the_read.size() ) ;

            AlignmentP frag_aln = init_alignment( size, size, 0, 0 ) ;

            std::string ref_for_mia = lifted ;
            for( size_t i = 0 ; i != ref_for_mia.length() ; ++i )
            {
                switch (toupper(ref_for_mia[i]))
                {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                        ref_for_mia[i] = toupper( ref_for_mia[i] ) ;
                        break ;
                    default:
                        ref_for_mia[i] = 'N' ;
                }
            }

            frag_aln->seq1 = ref_for_mia.c_str() ;
            frag_aln->len1 = ref_for_mia.size() ;
            frag_aln->seq2 = the_read.c_str() ;
            frag_aln->len2 = the_read.size() ;
            frag_aln->sg5 = 1 ;
            frag_aln->sg3 = 1 ;
            frag_aln->submat = submat ;
            pop_s1c_in_a( frag_aln ) ;
            pop_s2c_in_a( frag_aln ) ;
            dyn_prog( frag_aln ) ;

            pw_aln_frag pwaln ;
            max_sg_score( frag_aln ) ;			// ARGH!  This has a vital side-effect!!!
            find_align_begin( frag_aln ) ;  	//        And so has this...
            populate_pwaln_to_begin( frag_aln, &pwaln ) ;
            pwaln.start = frag_aln->abc;

            char *paln1 = aln_con, *paln2 = aln_ass ;
            int ass_pos = 0 ;
            while( ass_pos != (*s)->start && *paln1 && *paln2 )
            {
                if( *paln2 != '-' ) ass_pos++ ;
                ++paln1 ;
                ++paln2 ;
            }

            if( verbose >= 5 )
            {
                fprintf( stderr, "\n\naln.read: %s\naln.ref:  %s\nmatches:  ",
                         pwaln.frag_seq, pwaln.ref_seq ) ;

                const char* b = pwaln.ref_seq ;
                const char* pc = pwaln.frag_seq ;
                while( *b && *pc ) putc( *b++ == *pc++ ? '*' : ' ', stderr ) ;
                putc( '\n', stderr ) ;
                putc( '\n', stderr ) ;
            }

            cached_pwalns.push_back( cached_pwaln() ) ;
            cached_pwalns.back().start = pwaln.start ;
            cached_pwalns.back().ref_seq = pwaln.ref_seq ;
            cached_pwalns.back().frag_seq = pwaln.frag_seq ;

            std::string in_ref = lifted.substr( 0, pwaln.start ) ;
            in_ref.append( pwaln.ref_seq ) ;

            char *in_frag_v_ref = pwaln.frag_seq ;
            char *in_ass = maln->ref->seq + (*s)->start ;
            char *in_frag_v_ass = (*s)->seq ;

            if( verbose ) {
                if(*paln1!=in_ref[0]||*paln1=='-') fprintf( stderr, "huh? (R+%d) %.10s %.10s\n", pwaln.start, paln1, in_ref.c_str() ) ;
                if(*paln2!=in_ass[0]&&*paln2!='-') fprintf( stderr, "huh? (A+%d) %.10s %.10s\n", pwaln.start, paln2, in_ass ) ;
            }

            // iterate over alignment.  if we see something diagnosable
            // as contaminant, we mark that position as strong.
            while( ass_pos != (*s)->end +1 && *paln1 && *paln2 && !in_ref.empty() && *in_ass && *in_frag_v_ass && *in_frag_v_ref )
            {
                if( is_diagnostic( *paln1, *paln2 ) ) {
                    dp_list::iterator iter = l.find( ass_pos ) ;
                    if( iter == l.end() ) {
                        fprintf( stderr, "diagnostic site not found: %d\n", ass_pos ) ;
                    } else {
                        if( verbose >= 4 )
                            fprintf( stderr, "diagnostic pos.: %d %c(%c)/%c %c/%c",
                                    ass_pos, iter->second.consensus, in_ref[0], *in_frag_v_ref, *in_ass, *in_frag_v_ass ) ;
                        if( *in_frag_v_ref != *in_frag_v_ass )
                        {
                            if( verbose >= 4 ) fputs( " in disagreement.", stderr ) ;
                        } else {
                            bool maybe_clean = consistent( adna, iter->second.assembly, *in_frag_v_ass ) ;
                            bool maybe_dirt =  consistent( adna, iter->second.consensus,  *in_frag_v_ref ) ;

                            if( !maybe_clean && maybe_dirt && iter->second.strength == weak ) {
                                if( verbose >= 4 )
                                    fputs( " possible contaminant, upgraded to `effective'.", stderr ) ;
                                iter->second.contaminant = *in_frag_v_ref ;
                                iter->second.strength = effective ;
                            }
                        }
                    }
                    if( verbose >= 4 ) putc( '\n', stderr ) ;
                }

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
            if( verbose >= 4 ) fprintf( stderr, "\n" ) ;

            free_alignment( frag_aln ) ;
        }

        for( dp_list::iterator i = l.begin(), j = l.end() ; i != j ; )
        {
            dp_list::iterator k = i ;
            k++ ;
            if( i->second.strength == weak ) l.erase( i ) ;
            i=k ;
        }
        {
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
        }
        if( verbose >= 3 ) print_dp_list( stderr, l.begin(), l.end(), '\n' ) ;

        if( verbose >= 2 ) fputs( "Pass two: classifying fragments.\n", stderr ) ;
        std::deque< cached_pwaln >::const_iterator cpwaln = cached_pwalns.begin() ;
        for( const AlnSeqP *s = maln->AlnSeqArray ; s != maln->AlnSeqArray + maln->num_aln_seqs ; ++s, ++cpwaln )
        {
            whatsit klass = unknown ;
            whatsit klass2 = unknown ;
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
                    if( is_diagnostic( *paln1, *paln2 ) ) {
                        dp_list::const_iterator iter = l.find( ass_pos ) ;
                        if( iter != l.end() ) {
                            if( verbose >= 4 )
                                fprintf( stderr, "diagnostic pos. %s: %d %c(%c)/%c %c/%c",
                                        iter->second.strength == strong ? "(strong)" : "  (weak)",
                                        ass_pos, iter->second.consensus, in_ref[0], *in_frag_v_ref, *in_ass, *in_frag_v_ass ) ;
                            if( *in_frag_v_ref != *in_frag_v_ass )
                            {
                                if( verbose >= 4 ) fputs( " in disagreement.\n", stderr ) ;
                            }
                            else
                            {
                                bool maybe_clean = consistent( adna, iter->second.assembly, *in_frag_v_ass ) ;
                                bool maybe_dirt  = consistent( adna, iter->second.consensus,  *in_frag_v_ref ) ;

                                if( verbose >= 4 )
                                {
                                    fputs( maybe_dirt  ? " " : " in", stderr ) ;
                                    fputs( "consistent/", stderr ) ;
                                    fputs( maybe_clean ? "" : "in", stderr ) ;
                                    fputs( "consistent\n", stderr ) ;
                                }

                                update_class( klass2, votes2, maybe_clean, maybe_dirt && !maybe_clean ) ;
                                if( iter->second.strength == strong )
                                    update_class( klass, votes, maybe_clean, maybe_dirt ) ;
                            }
                        }
                    }

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

            Bfrags::const_iterator i = bfrags.find( (*s)->id ) ;
            Bfrags::const_iterator i2 = bfrags2.find( (*s)->id ) ;

            switch( (*s)->segment )
            {
                case 'b':
                    bfrags[ (*s)->id ] = std::make_pair( klass, votes ) ;
                    bfrags2[ (*s)->id ] = std::make_pair( klass2, votes2 ) ;
                    if( verbose >= 3 ) putc( '\n', stderr ) ;
                    break ;

                case 'f':
                    if( i == bfrags.end() )
                    {
                        fputs( (*s)->id, stderr ) ;
                        fputs( "/f is missing its back.\n", stderr ) ;
                    }
                    else
                    {
                        votes += i->second.second ;
                        klass = merge_whatsit( klass, i->second.first ) ;
                    }

                    if( i2 == bfrags2.end() )
                    {
                        fputs( (*s)->id, stderr ) ;
                        fputs( "/f is missing its back.\n", stderr ) ;
                    }
                    else
                    {
                        votes2 += i->second.second ;
                        klass2 = merge_whatsit( klass2, i->second.first ) ;
                    }

                case 'a':
                    if( verbose >= 2 ) fprintf( stderr, "%s is %s (%d votes)\n", (*s)->id, label[klass], votes ) ;
                    if( verbose >= 2 ) fprintf( stderr, "%s is %s (%d votes)\n", (*s)->id, label[klass2], votes2 ) ;
                    if( verbose >= 3 ) putc( '\n', stderr ) ;
                    summary[klass]++ ;
                    summary2[klass2]++ ;
                    break ;

                default:
                    fputs( "don't know how to handle fragment type ", stderr ) ;
                    putc( (*s)->segment, stderr ) ;
                    putc( '\n', stderr ) ;
            }
        }

        if( !mktable ) {
            int t = 0 ;
            for( dp_list::const_iterator i = l.begin(), e = l.end() ; i != e ; ++i )
                if( i->second.strength == strong ) t++ ;
            printf( "  strongly diagnostic positions: %d\n", t ) ;
        }
        print_results( summary, mktable ) ;
        if( !mktable ) printf( "  effectively diagnostic positions: %d\n", (int)l.size() ) ;
        else printf( "%d\t", (int)l.size() ) ;

        print_results( summary2, mktable ) ;
        putc( '\n', stdout ) ;

        free_map_alignment( maln ) ;
        free( aln_con ) ;
        free( aln_ass ) ;
    }
}

-}
