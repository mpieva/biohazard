module Bio.Glf (
        GlfSeq(..),
        GlfRec(..),
        enee_glf_file,
        enum_glf_file,
        enum_glf_handle
    ) where

import Bio.Iteratee
import Bio.Iteratee.Bgzf
import Control.Monad
import Data.Bits
import System.IO

import qualified Data.ByteString.Char8  as S
import qualified Data.Iteratee.ListLike as I


data GlfRec = SNP { glf_refbase :: {-# UNPACK #-} !Char
                  , glf_offset  :: {-# UNPACK #-} !Int
                  , glf_depth   :: {-# UNPACK #-} !Int
                  , glf_min_lk  :: {-# UNPACK #-} !Int
                  , glf_mapq    :: {-# UNPACK #-} !Int
                  , glf_lk      :: [Int] }
            | Indel { glf_refbase :: {-# UNPACK #-} !Char
                    , glf_offset  :: {-# UNPACK #-} !Int
                    , glf_depth   :: {-# UNPACK #-} !Int
                    , glf_min_lk  :: {-# UNPACK #-} !Int
                    , glf_mapq    :: {-# UNPACK #-} !Int
                    , glf_lk_hom1 :: {-# UNPACK #-} !Int
                    , glf_lk_hom2 :: {-# UNPACK #-} !Int
                    , glf_lk_het  :: {-# UNPACK #-} !Int
                    , glf_is_ins1 :: !Bool
                    , glf_is_ins2 :: !Bool
                    , glf_seq1    :: {-# UNPACK #-} !S.ByteString
                    , glf_seq2    :: {-# UNPACK #-} !S.ByteString }
    deriving Show

data GlfSeq = GlfSeq { glf_seqname :: {-# UNPACK #-} !S.ByteString
                     , glf_seqlen  :: {-# UNPACK #-} !Int }
    deriving Show


enee_glf_recs :: Monad m => Enumeratee S.ByteString [GlfRec] m b
enee_glf_recs = eneeCheckIfDone step
  where
    step  oit'       = I.isFinished >>= step' oit'

    step' oit'  True = return $ liftI oit'
    step' oit' False = do
        type_ref <- I.head
        let refbase = "XACMGRSVTWYHKDBN" !! fromIntegral (type_ref .&. 0xf)
        case type_ref `shiftR` 4 of
                0 -> return $ oit' $ EOF Nothing
                1 -> do r <- get_snp $ get_common (SNP refbase)
                        eneeCheckIfDone step . oit' $ Chunk [r]
                2 -> do r <- get_indel $ get_common (Indel refbase)
                        eneeCheckIfDone step . oit' $ Chunk [r]
                x -> fail $ "unknown GLF record #" ++ show x

    get_common f = return f
        `ap` (fromIntegral `liftM` endianRead4 LSB)
        `ap` (fromIntegral `liftM` endianRead3 LSB)
        `ap` (fromIntegral `liftM` I.head)
        `ap` (fromIntegral `liftM` I.head)

    get_snp f = f `ap` get_lk_arr
    get_lk_arr = replicateM 10 (fromIntegral `liftM` I.head)

    get_indel f = do
        f' <- f `ap` (fromIntegral `liftM` I.head)
                `ap` (fromIntegral `liftM` I.head)
                `ap` (fromIntegral `liftM` I.head)
        l1 <- getInt16le
        l2 <- getInt16le
        liftM2 (f' (l1 >= 0) (l2 >= 0)) (iGetString (abs l1)) (iGetString (abs l2))

    getInt16le = do i <- endianRead2 LSB
                    return $ if i > 0x7fff then fromIntegral i - 0x10000
                                           else fromIntegral i

enee_glf_seq :: Monad m => (GlfSeq -> Enumeratee [GlfRec] a m b) -> Enumeratee S.ByteString a m b
enee_glf_seq per_seq oit = do l <- endianRead4 LSB
                              s <- liftM2 GlfSeq (S.init `liftM` iGetString (fromIntegral l))
                                                 (fromIntegral `liftM` endianRead4 LSB)
                              enee_glf_recs ><> per_seq s $ oit

-- | Iterates over a GLF file.  In @get_glf_file per_seq per_file@, the
-- enumerator @per_file genome_name@, where @genome_name@ is the name
-- stored in the GLF header, is run once, then the enumeratee @per_seq
-- glfseq@ is iterated over the records in each sequence.
enee_glf_file :: Monad m => (GlfSeq -> Enumeratee [GlfRec] a m b)
                         -> (S.ByteString -> Enumerator a m b)
                         -> Enumeratee S.ByteString a m b
enee_glf_file per_seq per_file oit = do
    matched <- I.heads (S.pack "GLF\003")
    when (matched /= 4) (fail "GLF signature not found")
    nm <- endianRead4 LSB >>= iGetString . fromIntegral
    lift (per_file nm oit) >>= loop
  where
    -- loop :: Monad m => Iteratee a m b -> Iteratee S.ByteString m (Iteratee a m b)
    loop  it       = I.isFinished >>= loop' it
    loop' it  True = return it
    loop' it False = loop =<< enee_glf_seq per_seq it


-- | Enumerate the contents of a GLF file, apply suitable Enumeratees to
-- both sequences and records, resulting in an Enumerator of /whatever/,
-- typically output Strings or records...
--
-- This type is positively weird and I'm not entirely sure this is the
-- right way to go about it.
enum_glf_file :: MonadCatchIO m
              => FilePath
              -> (GlfSeq -> Enumeratee [GlfRec] a m b)
              -> (S.ByteString -> Enumerator a m b)
              -> Enumerator a m b
enum_glf_file fp per_seq per_file output =
    enumFile defaultBufSize fp >=> run $
    joinI $ decompressBgzf $
    enee_glf_file per_seq per_file output

enum_glf_handle :: MonadCatchIO m
                => Handle
                -> (GlfSeq -> Enumeratee [GlfRec] a m b)
                -> (S.ByteString -> Enumerator a m b)
                -> Enumerator a m b
enum_glf_handle hdl per_seq per_file output =
    enumHandle defaultBufSize hdl >=> run $
    joinI $ decompressBgzf $
    enee_glf_file per_seq per_file output

