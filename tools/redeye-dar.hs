{-# LANGUAGE FlexibleContexts #-}
-- Co-estimates aDNA damage with parameters for a simple genotype prior.
--
-- We want to estimate on only a subset of the genome.  For the time
-- being, this is by definition a subset of the large blocks of the
-- mappability track for the human genome (so this doesn't work for
-- other genomes).  To make this less crude, we need a differently
-- prepared input, but right now, input is one BAM file with read group
-- annotations.
--
-- We run the EM algorithm, repeatedly estimating one damage/error model
-- per read group and the global genotype parameters.  Convergence is
-- achieved when the changes in the damage models are sufficiently
-- small.  The damage/error model is one substitution matrix for each
-- position within a read near the ends, and one for what remains in the
-- middle.

import Bio.Adna
import Bio.Bam
import Bio.Bam.Pileup
import Bio.Genocall
import Bio.Genocall.Estimators          -- ( DivEst(..), tabulateSingle, estimateSingle, good_regions )
import Bio.Prelude
import Bio.Util.AD
import Data.Aeson.Encode.Pretty
import System.Console.GetOpt

import qualified Data.ByteString.Lazy.Char8     as L
import qualified Data.HashMap.Strict            as H
import qualified Data.Sequence                  as Z
import qualified Data.Vector                    as V
import qualified Data.Vector.Unboxed            as U
import qualified Data.Vector.Unboxed.Mutable    as M

data Conf = Conf {
    conf_output :: LazyBytes -> IO (),
    conf_report :: LazyBytes -> IO (),
    conf_params :: Parameters }

defaultConf :: Conf
defaultConf = Conf (L.hPutStrLn stdout) (\_ -> return ()) quietParameters

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o"  ["output"]     (ReqArg set_output "FILE") "Write output to FILE",
    Option "v"  ["verbose"]      (NoArg      set_verbose) "Print progress reports",
    Option "h?" ["help","usage"] (NoArg       disp_usage) "Print this message and exit" ]
    -- Missing here:  number of matrices?
  where
    set_verbose  c = return $ c { conf_report = L.hPutStrLn stderr, conf_params = debugParameters }
    set_output f c = return $ c { conf_output = L.writeFile f }

    disp_usage  _ = do pn <- getProgName
                       let blah = "Usage: " ++ pn ++ " [OPTION...] [LIBRARY-NAME...]"
                       putStrLn $ usageInfo blah options
                       exitSuccess

main :: IO ()
main = do
    (opts, files, errors) <- getOpt Permute options <$> getArgs
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Conf{..} <- foldl (>>=) (return defaultConf) opts

    -- For each iteration:  read the input, decompose, pileup.  The
    -- "prior" damage model is the usual 'SubstModel', the "posterior"
    -- damage model needs to be a mutable 'MSubstModel'.  We feed
    -- likelihoods into the div/het estimation (just as in
    -- 'redeye-pileup'), but also into a caller that will estimate
    -- damage.
    let iter sp0 mod0 = do ext_mod <- emIter sp0 mod0 files
                           conf_report $ encodePretty $ pop_separate ext_mod
                           if diffSubstMod mod0 (damage ext_mod) > 0.001
                               then iter (case point_est $ population ext_mod of
                                                [a,b] -> SinglePop a b) (damage ext_mod)
                               else return ext_mod

    final_model <- iter (SinglePop 0.001 0.002) (SubstModels H.empty)
    conf_output $ encodePretty (final_model :: ExtModel)

diffSubstMod :: SubstModels -> SubstModels -> Double
diffSubstMod (SubstModels m1) (SubstModels m2) =
    H.foldl' max 0 $
    H.map (either abs1 id) $
    H.unionWith diff1 (H.map Left m1) (H.map Left m2)
  where
    diff1 :: Either SubstModel a -> Either SubstModel a -> Either b Double
    diff1 (Left sm1) (Left sm2) = Right $ maximum $
            [ V.maximum $ V.zipWith diff2 (left_substs_fwd   sm1) (left_substs_fwd   sm2)
            ,                       diff2 (middle_substs_fwd sm1) (middle_substs_fwd sm2)
            , V.maximum $ V.zipWith diff2 (right_substs_fwd  sm1) (right_substs_fwd  sm2) ]

    diff2 :: Mat44D -> Mat44D -> Double
    diff2 (Mat44D u) (Mat44D v) = U.maximum $ U.map abs $ U.zipWith (-) u v

    abs1 :: SubstModel -> Double
    abs1 sm1 = V.maximum (V.map abs2 (left_substs_fwd   sm1)) `max`
                                abs2 (middle_substs_fwd sm1)  `max`
               V.maximum (V.map abs2 (right_substs_fwd  sm1))

    abs2 :: Mat44D -> Double
    abs2 (Mat44D v) = U.maximum v

-- One iteration of EM algorithm.  We go in with a substitution model
-- and het/div, we come out with new estimates for same.  We get het/div from
-- tabulation followed by numerical optimization.  For damage, we have to
-- compute posterior probabilities using the old model, then update the
-- damage matrices with pseudo counts.  (This amounts to a maximum
-- likelihood estimate for a weighted multinomial distribution.)
emIter :: SinglePop -> SubstModels -> [FilePath] -> IO ExtModel
emIter divest (SubstModels mod0) infiles =
        liftIO (mapM fresh_subst_model mod0 >>= newIORef)                         >>= \smodel ->
        concatInputs infiles >=> run                                                $ \hdr ->
        concatMapStreamM (decompose_dmg_from smodel)                               =$
        pileup                                                                     =$
        filterPilesWith (the_regions hdr)                                          =$
        mapStream  ( id &&& calls )                                                =$

        let div_estimation :: MonadIO m => Iteratee [(a, Calls)] m (DivEst,DivEst)
            div_estimation = mapStream snd                                         =$
                             tabulateSingle                                       >>=
                             liftIO . estimateSingle

            dmg_estimation :: MonadIO m => Iteratee [(Pile (Mat44D,MMat44D), Calls)] m (HashMap Bytes SubstModel)
            dmg_estimation = mapStreamM_ (\(p,c) ->
                                    liftIO . updateSubstModel p $
                                    single_pop_posterior divest
                                        (refix $ snp_refbase $ p_snp_pile c)
                                        (snp_gls $ p_snp_pile c))                  >>
                             liftIO (readIORef smodel)                            >>=
                             liftIO . mapM (freezeSubstModel . snd)

        in (\((d1,d2),m) -> ExtModel d1 (Just d2) (SubstModels m))
           <$> zipStreams div_estimation dmg_estimation
  where
    the_regions hdr = sort [ Region (Refseq $ fromIntegral ri) p (p+l)
                           | (ch, p, l) <- good_regions
                           , let Just ri = Z.findIndexL ((==) ch . sq_name) (meta_refs hdr) ]

    refix ref = U.fromListN 16 [0,0,2,0,5,0,0,0,9,0,0,0,0,0,0,0] U.! fromIntegral (unNs ref)


filterPilesWith :: Monad m => [Region] -> Enumeratee [Pile a] [Pile a] m b
filterPilesWith = unfoldConvStream go
  where
    go [    ] = ([],[]) <$ skipToEof
    go (r:rs) = do mp <- peekStream
                   case mp of
                        Just p | (p_refseq p, p_pos p) <  (refseq r, start r) -> headStream >> go (r:rs)
                               | (p_refseq p, p_pos p) >= (refseq r, end   r) -> go rs
                               | otherwise                                    -> (\x -> (r:rs, [x])) <$> headStream
                        Nothing                                               -> return ([],[])


-- Probabilistically count substitutions.  We infer from posterior
-- genotype probabilities what the base must have been, then count
-- substitutions from that to the actual base.
updateSubstModel :: Pile ( Mat44D, MMat44D ) -> U.Vector Prob -> IO ()
updateSubstModel pile postp = mapM_ count_base bases
  where
    -- Posterior probalities of the haploid base before damage
    -- @P(H) = \sum_{G} P(H|G) P(G|D)@
    pH_A = fromProb $ postp U.! 0 + 0.5 * ( postp U.! 1 + postp U.! 3 + postp U.! 6 )
    pH_C = fromProb $ postp U.! 2 + 0.5 * ( postp U.! 1 + postp U.! 4 + postp U.! 7 )
    pH_G = fromProb $ postp U.! 5 + 0.5 * ( postp U.! 3 + postp U.! 4 + postp U.! 8 )
    pH_T = fromProb $ postp U.! 9 + 0.5 * ( postp U.! 6 + postp U.! 7 + postp U.! 8 )

    -- XXX This ignores map quality.  I don't think it matters.
    bases = map snd $ uncurry (++) $ p_snp_pile pile

    -- P(H:->X) = P(H|X)
    --          = P(X|H) P(H) / P(X)
    --          = P(X|H) P(H) / \sum_H' P(X|H') P(H')
    --
    -- We get P(X|H) from the old substitution model.
    -- Fortunately, it's actually available.

    count_base b = do nudge (snd $ db_dmg b) (nucA :-> db_call b) (pHX_A/pHX)
                      nudge (snd $ db_dmg b) (nucC :-> db_call b) (pHX_C/pHX)
                      nudge (snd $ db_dmg b) (nucG :-> db_call b) (pHX_G/pHX)
                      nudge (snd $ db_dmg b) (nucT :-> db_call b) (pHX_T/pHX)
      where
        pHX = pHX_A + pHX_C + pHX_G + pHX_T
        pHX_A = (fst (db_dmg b) `bang` nucA :-> db_call b) * pH_A
        pHX_C = (fst (db_dmg b) `bang` nucC :-> db_call b) * pH_C
        pHX_G = (fst (db_dmg b) `bang` nucG :-> db_call b) * pH_G
        pHX_T = (fst (db_dmg b) `bang` nucT :-> db_call b) * pH_T



calls :: Pile (Mat44D,a) -> Calls
calls pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !s = simple_snp_call   $ map (second (fmap fst)) $ uncurry (++) $ p_snp_pile pile
    !i = simple_indel_call $ map (second (second (map (fmap fst)))) $ p_indel_pile pile
    --
    -- XXX this should be a cmdline option, if we ever look at qualities again
    -- fq = min 1 . (*) 1.333 . fromQual
    -- fq = fromQual


{-# INLINE decompose_dmg_from #-}
decompose_dmg_from :: IORef (HashMap Bytes (SubstModel, MSubstModel)) -> BamRaw -> IO [ PosPrimChunks (Mat44D, MMat44D) ]
decompose_dmg_from ref raw = do
    hm <- readIORef ref
    let rg = extAsString "RG" (unpackBam raw)
    model <- case H.lookup rg hm of
                Just mm -> return mm
                Nothing -> do mm <- fresh_subst_model model0
                              writeIORef ref $! H.insert rg mm hm
                              return mm
    return $ decompose (from_mmodel model) raw

  where
    from_mmodel (m,mm) i r
        | i >= 0 &&   i  <  V.length  (left_substs_fwd m) && not r
                = ( V.unsafeIndex (left_substs_fwd   m)   i,    V.unsafeIndex (left_substs_fwd  mm)   i )

        | i <  0 && (-i) <= V.length (right_substs_fwd m) && not r
                = ( V.unsafeIndex (right_substs_fwd  m) (-i-1), V.unsafeIndex (right_substs_fwd mm) (-i-1) )

        | not r = ( middle_substs_fwd m, middle_substs_fwd mm )

        | i >= 0 &&   i  <  V.length  (left_substs_rev m)
                = ( V.unsafeIndex (left_substs_rev   m)   i,    V.unsafeIndex (left_substs_rev  mm)   i )

        | i <  0 && (-i) <= V.length (right_substs_rev m)
                = ( V.unsafeIndex (right_substs_rev  m) (-i-1), V.unsafeIndex (right_substs_rev mm) (-i-1) )

        | True  = ( middle_substs_rev m, middle_substs_rev mm )


    initmat = Mat44D $ U.fromListN 16 [ 0.91, 0.03, 0.03, 0.03
                                      , 0.03, 0.91, 0.03, 0.03
                                      , 0.03, 0.03, 0.91, 0.03
                                      , 0.03, 0.03, 0.03, 0.91 ]

    model0 = SubstModel { left_substs_fwd   = V.replicate 12 initmat
                        , middle_substs_fwd =                initmat
                        , right_substs_fwd  = V.replicate 12 initmat
                        , left_substs_rev   = V.replicate 12 initmat
                        , middle_substs_rev =                initmat
                        , right_substs_rev  = V.replicate 12 initmat }

fresh_subst_model :: SubstModel -> IO (SubstModel, MSubstModel)
fresh_subst_model m = (,) m <$> m'
  where
    m' = SubstModel <$> V.mapM nullmat (left_substs_fwd   m)
                    <*>        nullmat (middle_substs_fwd m)
                    <*> V.mapM nullmat (right_substs_fwd  m)
                    <*> V.mapM nullmat (left_substs_rev   m)
                    <*>        nullmat (middle_substs_rev m)
                    <*> V.mapM nullmat (right_substs_rev  m)

    nullmat = const $ MMat44D <$> M.replicate 16 (0::Double)

