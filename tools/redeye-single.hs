-- Genotype calling for a single individual.  The only parameters needed
-- are the (prior) probabilities for a heterozygous or homozygous variant.
--
-- (Could this be extended easily into a caller for a homogenous
-- population?  Individuals would be assumed to be from a homogenous
-- population and not closely related (i.e. not family).  In this case,
-- the prior is the allele frequency spectrum, the call would be the
-- set(!) of genotypes that has maximum posterior probability.  It's not
-- quite clear how to compute it efficiently, though.)
--
-- What's the output format?  Fasta could be useful in limited
-- circumstances, VCF or BCF might be useful, too.  Or maybe VCF
-- restricted to variant sites.  Or VCF restricted to sites not known to
-- be reference.  People will come up with filters for sure...
--
-- Sometimes we might want to fuse this with simple calling
-- for one individual, in which case the output becomes VCF(?), BCF or
-- maybe Fasta/Fastq.  People will want to filter, too.  When fused, we
-- might also want to estimate parameters.  And when not fused, we read
-- the Avro container back in.
--
-- Maybe do it separately first?

data Conf = Conf {
    conf_output      :: Maybe Output,
    conf_sample      :: S.ByteString,
    conf_ploidy      :: S.ByteString -> Int,
    conf_report      :: String -> IO (),
    conf_prior_het   :: Prob Double,
    conf_prior_indel :: Prob Double }

defaultConf :: Conf
defaultConf = Conf Nothing "John_Doe" (const 2) (\_ -> return ())
                   (qualToProb $ Q 30) (qualToProb $ Q 45)


type OIter = Conf -> Refs -> Iteratee [Calls] IO ()

output_fasta :: FilePath -> (OIter -> IO r) -> IO r
output_fasta fn k = if fn == "-" then k (fa_out stdout)
                                 else withFile fn WriteMode $ k . fa_out
  where
    fa_out :: Handle -> Conf -> Refs -> Iteratee [Calls] IO ()
    fa_out hdl Conf{..} refs =
            by_groups p_refseq (\rs out -> do
                    let sname = sq_name $ getRef refs rs
                    out' <- lift $ enumPure1Chunk [S.concat [">", conf_sample, "--", sname]] out
                    convStream (do callz <- headStream
                                   let s1 = format_snp_call conf_prior_het callz
                                   S.append s1 <$> format_indel_call conf_prior_indel callz)
                          =$ collect_lines out') =$
            mapStreamM_ (S.hPut hdl . (flip S.snoc '\n'))

-- | Formatting a SNP call.  If this was a haplopid call (four GL
-- values), we pick the most likely base and pass it on.  If it was
-- diploid, we pick the most likely dinucleotide and pass it on.

format_snp_call :: Prob Double -> Calls -> S.ByteString
format_snp_call p cs
    | U.length gl ==  4 = S.take 1 $ S.drop (maxQualIndex gl) hapbases
    | U.length gl == 10 = S.take 1 $ S.drop (maxQualIndex $ U.zipWith (*) ps gl) dipbases
    | otherwise = error "Thou shalt not try to format_snp_call unless thou madeth a haploid or diploid call!"
  where
    gl = p_snp_pile cs
    ps = U.fromListN 10 [p,1,p,1,1,p,1,1,1,p]
    dipbases = "NAMCRSGWYKT"
    hapbases = "NACGT"

-- | Formatting an Indel call.  We pick the most likely variant and
-- pass its sequence on.  Then we drop incoming calls that should be
-- deleted according to the chosen variant.  Note that this will blow up
-- unless the call was done assuming a haploid genome (which is
-- guaranteeed /in this program/)!

format_indel_call :: Monad m => Prob Double -> Calls -> Iteratee [Calls] m S.ByteString
format_indel_call p cs
    | U.length gl0 == nv                  = go gl0
    | U.length gl0 == nv * (nv+1) `div` 2 = go homs
    | otherwise = error "Thou shalt not try to format_indel_call unless thou madeth a haploid or diploid call!"
  where
    (gl0,vars) = p_indel_pile cs
    !nv   = length vars
    !homs = U.fromListN nv [ gl0 U.! (i*(i+1) `div` 2 -1) | i <- [1..nv] ]

    go gl = I.dropWhile skip >> return (S.pack $ show $ U.toList ins)
      where
        eff_gl = U.fromList $ zipWith adjust (U.toList gl) vars
        adjust q (IndelVariant ds (V_Nuc is)) = if ds == 0 && U.null is then q else p * q

        IndelVariant del (V_Nuc ins) = ( IndelVariant 0 (V_Nuc U.empty) : vars ) !! maxQualIndex eff_gl
        skip ocs  = p_refseq ocs == p_refseq cs && p_pos ocs < p_pos cs + del

maxQualIndex :: U.Vector (Prob Double) -> Int
maxQualIndex vec = case U.ifoldl' step (0, 0, 0) vec of
    (!i, !m, !m2) -> if m / m2 > 2 then i else 0
  where
    step (!i,!m,!m2) j v = if v >= m then (j+1,v,m) else (i,m,m2)

collect_lines :: Monad m => Enumeratee S.ByteString [S.ByteString] m r
collect_lines = eneeCheckIfDone (liftI . go S.empty)
  where
    go acc k (EOF  mx) = idone (k $ Chunk [acc]) $ EOF mx
    go acc k (Chunk s) = case S.splitAt 60 (acc `S.append` s) of
                            (left, right) | S.null right -> liftI $ go left k
                                          | otherwise    -> eneeCheckIfDone (liftI . go right) . k $ Chunk [left]

