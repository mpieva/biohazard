{-# LANGUAGE FlexibleContexts #-}
-- Estimates aDNA damage.  Crude first version.
--
-- - Read or subsample a BAM file, make compact representation of the reads.
-- - Compute likelihood of each read under simple model of
--   damage, error/divergence, contamination.
--
-- For the fitting, we simplify radically: ignore sequencing error,
-- assume damage and simple, symmetric substitutions which subsume error
-- and divergence.
--
-- Trying to compute symbolically is too much, the high power terms get
-- out of hand quickly, and we get mixed powers of \lambda and \kappa.
-- The fastest version so far uses the cheap implementation of automatic
-- differentiation in AD.hs together with the Hager-Zhang method from
-- package nonlinear-optimization.  BFGS from hmatrix-gsl takes longer
-- to converge.  Didn't try an actual Newton iteration (yet?), AD from
-- package ad appears slower.
--
-- If I include parameters, whose true value is zero, the transformation
-- to the log-odds-ratio doesn't work, because then the maximum doesn't
-- exist anymore.  For many parameters, zero makes sense, but one
-- doesn't.  A different transformation ('sigmoid2'/'isigmoid2'
-- below) allows for an actual zero (but not an actual one), while
-- avoiding ugly boundary conditions.  That appears to work well.
--
-- The current hack assumes all molecules have an overhang at both ends,
-- then each base gets deaminated with a position dependent probability
-- following a geometric distribution.  If we try to model a fraction of
-- undeaminated molecules (a contaminant) in addition, this fails.  To
-- rescue the idea, I guess we must really decide if the molecule has an
-- overhang at all (probability 1/2) at each end, then deaminate it.

import Bio.Adna
import Bio.Bam
import Bio.Bam.Pileup
import Bio.Genocall
import Bio.Genocall.Estimators          ( DivEst(..), tabulateSingle, estimateSingle, good_regions )
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
import qualified Text.PrettyPrint.Leijen.Text   as P

data Conf = Conf {
    conf_report :: String -> IO (),
    conf_params :: Parameters }

defaultConf :: Conf
defaultConf = Conf (\_ -> return ()) quietParameters

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "v"  ["verbose"]      (NoArg      set_verbose) "Print progress reports",
    Option "h?" ["help","usage"] (NoArg       disp_usage) "Print this message and exit" ]
  where
    set_verbose  c = return $ c { conf_report = hPutStrLn stderr, conf_params = debugParameters }

    disp_usage  _ = do pn <- getProgName
                       let blah = "Usage: " ++ pn ++ " [OPTION...] [LIBRARY-NAME...]"
                       putStrLn $ usageInfo blah options
                       exitSuccess

main :: IO ()
main = do
    (opts, files, errors) <- getOpt Permute options <$> getArgs
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Conf{..} <- foldl (>>=) (return defaultConf) opts

    -- Meh.  New concept.  We'll operate repeatedly on a small set of
    -- regions.  For the time being, that region will have to be put
    -- into a suitable file, so here we expect precisely one file.  Read
    -- group annotations have to be present, we will apply one damage
    -- model per read group.  Have to decide when we're done.

    -- For each iteration:  read the input, decompose, pileup.  The
    -- "prior" damage model is the usual 'SubstModel', the "posterior"
    -- damage model needs to be a mutable 'MSubstModel'.  We feed
    -- likelihoods into the div/het estimation (just as in
    -- 'redeye-pileup'), but also into a caller that will estimate
    -- damage.
    let iter sp0 mod0 = do ((u,v), model1) <- emIter sp0 mod0 conf_report files
                           L.putStrLn $ encodePretty (u,v)
                           -- pprint $ H.toList model1
                           iter (case point_est u of [d,h] -> SinglePop d h) model1

    iter (SinglePop 0.001 0.002) H.empty


-- One iteration of EM algorithm.  We go in with a substitution model
-- and het/div, we come out with new estimates.  We get het/div from
-- tabulation followed by estimation, as before.  For damage, we have to
-- compute posterior probabilities using the old model, then update the
-- damage probabilistically.
emIter :: SinglePop -> HashMap Bytes SubstModel -> (String -> IO ())
       -> [FilePath] -> IO ((DivEst, DivEst), HashMap Bytes SubstModel)
emIter divest mod0 report infiles =
        liftIO (mapM fresh_subst_model mod0 >>= newIORef)                         >>= \smodel ->
        concatInputs infiles >=> run                                                $ \hdr ->
        concatMapStreamM (decompose_dmg_from smodel)                               =$
        progressPos (\(rs,p,_)->(rs,p)) "GT call at" (meta_refs hdr) 0x8000 report =$
        pileup                                                                     =$
        filterPilesWith (the_regions hdr)                                          =$
        mapStream  ( id &&& calls )                                                =$
        zipStreams ( mapStream snd                                      =$
                     tabulateSingle                                     >>=
                     liftIO . estimateSingle )
                   ( mapStreamM_ (uncurry (updateSubstModel divest))    >>
                     liftIO (readIORef smodel)                          >>=
                     liftIO . mapM (freezeSubstModel . snd) )
  where
    the_regions hdr = sort [ Region (Refseq $ fromIntegral ri) p (p+l)
                           | (ch, p, l) <- good_regions
                           , let Just ri = Z.findIndexL ((==) ch . sq_name) (meta_refs hdr) ]


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
updateSubstModel :: SinglePop -> Pile ( Mat44D, MMat44D ) -> Calls -> IO ()
updateSubstModel divest pile cs = mapM_ count_base bases
  where
    postp = single_pop_posterior divest (snp_refbase (p_snp_pile cs)) (snp_gls (p_snp_pile cs))

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



freezeSubstModel :: MSubstModel -> IO SubstModel
freezeSubstModel mm = do
    new_left   <- V.zipWithM freezeMats (left_substs_fwd   mm) (right_substs_rev  mm)
    new_middle <-            freezeMats (middle_substs_fwd mm) (middle_substs_rev mm)
    new_right  <- V.zipWithM freezeMats (right_substs_fwd  mm) (left_substs_rev   mm)

    return $ SubstModel new_left new_middle new_right
                        ( V.map compl_mat new_left   )
                              ( compl_mat new_middle )
                        ( V.map compl_mat new_right  )
  where
    -- We take two matrices (forward and reverse strand model), flip the
    -- reverse one around, normalize for the total ocunt and freeze the
    -- result into a substitution matrix.  The two matrices must
    -- correspond, so the rev-strand vector has to be reversed.
    freezeMats :: MMat44D -> MMat44D -> IO Mat44D
    freezeMats (MMat44D vv) (MMat44D ww) = do
        v <-             Mat44D <$> U.freeze vv
        w <- compl_mat . Mat44D <$> U.freeze ww
        return . Mat44D $ U.fromListN 16
                [ ((v `bang` x :-> y) + (w `bang` x :-> y)) / s
                | x <- range (nucA, nucT)
                , let s = sum [ (v `bang` x :-> y) + (w `bang` x :-> y)
                              | y <- range (nucA, nucT) ]
                , y <- range (nucA, nucT) ]

    -- Huh, I thought this should be transposed.  But *that* doesn't work.
    compl_mat :: Mat44D -> Mat44D
    compl_mat v = Mat44D $ U.fromListN 16
                    [ v `bang` compl x :-> compl y
                    | y <- range (nucA, nucT)
                    , x <- range (nucA, nucT) ]

calls :: Pile ( Mat44D, b ) -> Calls
calls pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !s = simple_snp_call   $ map (second (fmap fst)) $ uncurry (++) $ p_snp_pile pile
    !i = simple_indel_call $ map (second (second (map (fmap fst)))) $ p_indel_pile pile


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

