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
import Bio.Genocall.Estimators
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
    conf_params :: Parameters,
    conf_length :: Int,
    conf_eps    :: Double }

defaultConf :: Conf
defaultConf = Conf (L.hPutStrLn stdout) (\_ -> return ()) quietParameters 16 1.0E-6

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o"  ["output"]     (ReqArg set_output "FILE") "Write output to FILE (stdout)",
    Option "l"  ["model-length"] (ReqArg   set_len "NUM") "Set size of subst. model to NUM (16)",
    Option "e"  ["precision"]    (ReqArg  set_prec "NUM") "Set precision for fit to NUM (1E-6)",
    Option "v"  ["verbose"]      (NoArg      set_verbose) "Print progress reports",
    Option "h?" ["help","usage"] (NoArg       disp_usage) "Print this message and exit" ]
  where
    set_verbose  c = return $ c { conf_report = L.hPutStrLn stderr, conf_params = debugParameters }
    set_output f c =                    return $ c { conf_output = L.writeFile f }
    set_len    a c = readIO a >>= \x -> return $ c { conf_length = x }
    set_prec   a c = readIO a >>= \x -> return $ c { conf_eps    = x }

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
    let iter sp0 mod0 = do (((de1,de2),mod1),syms) <- emIter conf_length sp0 mod0 files
                           conf_report $ encodePretty de2
                           if diffSubstMod mod0 mod1 > conf_eps
                               then iter (case point_est de1 of [a,b] -> SinglePop a b) mod1
                               else return . ExtModel de1 (Just de2) . SubstModels
                                           $ H.map ((mod1 V.!) . fromDmgToken) syms

    final_model <- iter (SinglePop 0.001 0.002) V.empty
    conf_output $ encodePretty (final_model :: ExtModel)

diffSubstMod :: V.Vector SubstModel -> V.Vector SubstModel -> Double
diffSubstMod v1 v2 =
    V.foldl' max 0 (V.zipWith diff1 v1 v2) `max`
    V.foldl' max 0 (V.map abs1 (V.drop (V.length v2) v1)) `max`
    V.foldl' max 0 (V.map abs1 (V.drop (V.length v1) v2))
  where
    diff1 :: SubstModel -> SubstModel -> Double
    diff1 sm1 sm2 = maximum $
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
emIter :: Int -> SinglePop -> V.Vector SubstModel -> [FilePath] -> IO (((DivEst,DivEst), V.Vector SubstModel), HashMap Bytes DmgToken)
emIter msize divest mod0 infiles =
        liftIO (newIORef V.empty)                                                 >>= \mmod ->
        liftIO (newIORef H.empty)                                                 >>= \symtab ->
        concatInputs infiles >=> run                                                $ \hdr ->
        concatMapStreamM (decompose_dmg_from symtab)                               =$
        pileup                                                                     =$
        filterPilesWith (the_regions hdr)                                          =$
        mapStream ( id &&& calls mod0 )                                            =$

        let div_estimation :: MonadIO m => Iteratee [(a, Calls)] m (DivEst,DivEst)
            div_estimation = mapStream snd                                         =$
                             tabulateSingle                                       >>=
                             liftIO . estimateSingle

            dmg_estimation :: MonadIO m => Iteratee [(Pile, Calls)] m (V.Vector SubstModel)
            dmg_estimation = mapStreamM_ (\(p,c) ->
                                    liftIO . updateSubstModel msize mod0 mmod p $
                                    single_pop_posterior divest
                                        (refix $ snp_refbase $ p_snp_pile c)
                                        (snp_gls $ p_snp_pile c))                  >>
                             liftIO (readIORef mmod)                              >>=
                             liftIO . V.mapM freezeSubstModel

        in (,) <$> zipStreams div_estimation dmg_estimation
               <*> liftIO (readIORef symtab)
  where
    the_regions hdr = sort [ Region (Refseq $ fromIntegral ri) p (p+l)
                           | (ch, p, l) <- good_regions
                           , let Just ri = Z.findIndexL ((==) ch . sq_name) (meta_refs hdr) ]

    refix ref = U.fromListN 16 [0,0,2,0,5,0,0,0,9,0,0,0,0,0,0,0] U.! fromIntegral (unNs ref)


filterPilesWith :: Monad m => [Region] -> Enumeratee [Pile] [Pile] m b
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

updateSubstModel :: Int -> V.Vector SubstModel -> IORef (V.Vector MSubstModel) -> Pile -> U.Vector Prob -> IO ()
updateSubstModel msize mods0 vmods1 pile postp = case p_snp_pile pile of
    -- XXX This ignores map quality.  I don't think it matters.
    (basesF, basesR) -> do mapM_ (count_base False . snd) basesF
                           mapM_ (count_base  True . snd) basesR
  where
    -- Posterior probalities of the haploid base before damage
    -- @P(H) = \sum_{G} P(H|G) P(G|D)@
    pH_A = fromProb $ postp U.! 0 + 0.5 * ( postp U.! 1 + postp U.! 3 + postp U.! 6 )
    pH_C = fromProb $ postp U.! 2 + 0.5 * ( postp U.! 1 + postp U.! 4 + postp U.! 7 )
    pH_G = fromProb $ postp U.! 5 + 0.5 * ( postp U.! 3 + postp U.! 4 + postp U.! 8 )
    pH_T = fromProb $ postp U.! 9 + 0.5 * ( postp U.! 6 + postp U.! 7 + postp U.! 8 )

    -- bases = map snd $ uncurry (++) $ p_snp_pile pile

    -- P(H:->X) = P(H|X)
    --          = P(X|H) P(H) / P(X)
    --          = P(X|H) P(H) / \sum_H' P(X|H') P(H')
    --
    -- We get P(X|H) from the old substitution model.
    -- Fortunately, it's actually available.

    count_base str b = do
        let old_mat = case mods0 V.!? fromDmgToken (db_dmg_tk b) of
                        Nothing -> initmat
                        Just sm -> lookupSubstModel sm (db_dmg_pos b) str

        new_mat <- do mods1 <- readIORef vmods1
                      sm <- case mods1 V.!? fromDmgToken (db_dmg_tk b) of
                            Nothing -> do m <- new_mmodel
                                          writeIORef vmods1 (V.snoc mods1 m)
                                          return m
                            Just  m -> return m
                      return $ lookupSubstModel sm (db_dmg_pos b) str

        let pHX = pHX_A + pHX_C + pHX_G + pHX_T
            pHX_A = (old_mat `bang` nucA :-> db_call b) * pH_A
            pHX_C = (old_mat `bang` nucC :-> db_call b) * pH_C
            pHX_G = (old_mat `bang` nucG :-> db_call b) * pH_G
            pHX_T = (old_mat `bang` nucT :-> db_call b) * pH_T

        nudge new_mat (nucA :-> db_call b) (pHX_A/pHX)
        nudge new_mat (nucC :-> db_call b) (pHX_C/pHX)
        nudge new_mat (nucG :-> db_call b) (pHX_G/pHX)
        nudge new_mat (nucT :-> db_call b) (pHX_T/pHX)


    new_mmodel = SubstModel <$> V.replicateM msize nullmat
                            <*>                    nullmat
                            <*> V.replicateM msize nullmat
                            <*> V.replicateM msize nullmat
                            <*>                    nullmat
                            <*> V.replicateM msize nullmat

    nullmat = MMat44D <$> M.replicate 16 (0::Double)

calls :: V.Vector SubstModel -> Pile -> Calls
calls dmg pile = pile { p_snp_pile = s, p_indel_pile = i }
  where
    !s = simple_snp_call   get_dmg $ p_snp_pile pile
    !i = simple_indel_call get_dmg $ p_indel_pile pile

    get_dmg :: DmgToken -> Int -> Bool -> Mat44D
    get_dmg (DmgToken dt) di ds = case dmg V.!? dt of
        Nothing -> initmat                      -- not found, happens in first round
        Just sm -> lookupSubstModel sm di ds

initmat :: Mat44D
initmat = Mat44D $ U.fromListN 16 [ 0.91, 0.03, 0.03, 0.03
                                  , 0.03, 0.91, 0.03, 0.03
                                  , 0.03, 0.03, 0.91, 0.03
                                  , 0.03, 0.03, 0.03, 0.91 ]

{-# INLINE decompose_dmg_from #-}
decompose_dmg_from :: IORef (HashMap Bytes DmgToken) -> BamRaw -> IO [PosPrimChunks]
decompose_dmg_from ref raw = do
    hm <- readIORef ref
    let rg = extAsString "RG" (unpackBam raw)
    token <- case H.lookup rg hm of
                Just tk -> return tk
                Nothing -> do let tk = DmgToken $ H.size hm
                              writeIORef ref $! H.insert rg tk hm
                              return tk
    return $ decompose token raw


