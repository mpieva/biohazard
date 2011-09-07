module Bio.File.Glf (
    GlfSeq(..),
    GlfRec(..),
    enee_glf_file,
    enum_glf_file,
    enum_glf_handle
                    ) where

import Bio.File.Bgzf
import Bio.Util
import Control.Monad
import Control.Monad.CatchIO
import Control.Monad.Trans.Class
import Data.Bits
import System.IO

import qualified Data.ByteString.Char8  as S
import qualified Data.Iteratee.Binary   as I
import qualified Data.Iteratee.IO       as I
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


enee_glf_recs :: Monad m => I.Enumeratee [GlfRec] a m b -> I.Enumeratee S.ByteString a m b
enee_glf_recs enee oit = step (enee oit)
  where
    -- step :: Monad m => (I.Iteratee [GlfRec] m (I.Iteratee a m b)) -> I.Iteratee S.ByteString m (I.Iteratee a m b)
    step  oit'       = I.isFinished >>= step' oit'

    -- step' :: Monad m => (I.Iteratee [GlfRec] m (I.Iteratee a m b)) -> Bool -> I.Iteratee S.ByteString m (I.Iteratee a m b)
    step' oit'  True = lift $ I.run =<< I.enumEof oit'
    step' oit' False = do
        type_ref <- I.head
        let refbase = "XACMGRSVTWYHKDBN" !! fromIntegral (type_ref .&. 0xf)
        case type_ref `shiftR` 4 of
                0 -> lift $ I.run =<< I.enumEof oit'
                1 -> do r <- get_snp $ get_common (SNP refbase)
                        lift (I.enumPure1Chunk [r] oit') >>= step
                2 -> do r <- get_indel $ get_common (Indel refbase)
                        lift (I.enumPure1Chunk [r] oit') >>= step
                x -> fail $ "unknown GLF record #" ++ show x

    get_common f = return f
        `ap` (fromIntegral `liftM` I.endianRead4 I.LSB)
        `ap` (fromIntegral `liftM` I.endianRead3 I.LSB)
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
        liftM2 (f' (l1 >= 0) (l2 >= 0)) (getString (abs l1)) (getString (abs l2))

    getInt16le = do i <- I.endianRead2 I.LSB
                    return $ if i > 0x7fff then fromIntegral i - 0x10000
                                           else fromIntegral i

enee_glf_seq :: Monad m => (GlfSeq -> I.Enumeratee [GlfRec] a m b) -> I.Enumeratee S.ByteString a m b
enee_glf_seq per_seq oit = do l <- I.endianRead4 I.LSB
                              s <- liftM2 GlfSeq (S.init `liftM` getString (fromIntegral l))
                                                 (fromIntegral `liftM` I.endianRead4 I.LSB)
                              enee_glf_recs (per_seq s) oit

-- | Iterates over a GLF file.  In @get_glf_file per_seq per_file@, the
-- enumerator @per_file genome_name@, where @genome_name@ is the name
-- stored in the GLF header, is run once, then the enumeratee @per_seq
-- glfseq@ is iterated over the records in each sequence.
enee_glf_file :: Monad m => (GlfSeq -> I.Enumeratee [GlfRec] a m b)
                         -> (S.ByteString -> I.Enumerator a m b)
                         -> I.Enumeratee S.ByteString a m b
enee_glf_file per_seq per_file oit = do 
    matched <- I.heads (S.pack "GLF\003")
    when (matched /= 4) (fail "GLF signature not found")
    nm <- I.endianRead4 I.LSB >>= getString . fromIntegral
    lift (per_file nm oit) >>= loop
  where
    -- loop :: Monad m => I.Iteratee a m b -> I.Iteratee S.ByteString m (I.Iteratee a m b)
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
              -> (GlfSeq -> I.Enumeratee [GlfRec] a m b)
              -> (S.ByteString -> I.Enumerator a m b)
              -> I.Enumerator a m b
enum_glf_file fp per_seq per_file output =
    I.enumFile I.defaultBufSize fp >=> I.run $
    I.joinI $ decompress $
    enee_glf_file per_seq per_file output

enum_glf_handle :: MonadCatchIO m
                => Handle 
                -> (GlfSeq -> I.Enumeratee [GlfRec] a m b) 
                -> (S.ByteString -> I.Enumerator a m b)
                -> I.Enumerator a m b
enum_glf_handle hdl per_seq per_file output =
    I.enumHandle I.defaultBufSize hdl >=> I.run $
    I.joinI $ decompress $
    enee_glf_file per_seq per_file output

