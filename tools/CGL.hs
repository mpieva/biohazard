module CGL where

-- "Compact Genotype Likelihoods"
--
-- Serialize the results from genotype calling in a sensible way.  A
-- file consists of 'Block's, each 'Block' comes with a header
-- identifying the reference sequence, and many 'Variant' records.
-- The file itself gets a minimal header to facilitate versioning, but
-- not much else.  We might want to compress 'Block's.
--
-- The original idea is to create an AVRO container file, but that
-- requires another layer of blocking, which I'm not going to implement
-- right now.  The encoding should be nearly valid AVRO, though.

encodeCGL :: Monad m => S.ByteString -> Refs -> Enumeratee [Calls] S.ByteString m a
encodeCGL key refs = enumBuilder (byteString "CGL\0" <> aString key) >=> oneBlock
  where
    oneBlock out = do
        mc1 <- peek
        case mc1 of
            Nothing -> return out
            Just c1 -> tailBlock out (p_refseq c1) (p_pos c1) (p_pos c1) 1024 []

    tailBlock out !rs !p0 !po !n acc out = do
        mc <- peek
        case mc of
            Just c | rs == p_refseq c1 && po+1 == p_pos c1 && n > 0 ->
                tailBlock out rs p0 (po+1) (n-1) $ (p_snp c, p_indel c) : acc

            _ -> let bld = mconcat $
                        aString (sq_name (getRef refs (p_refseq c1))) :
                        anInt (p_pos c1) :
                        anInt (length acc) :
                        map twoCalls (reverse acc)
                 in enumBuilder bld out >>= oneBlock

enumBuilder :: Monad m => Builder -> Enumeratee S.ByteString m a
enumBuilder = enumList . toChunks . toLazyByteString

anInt -- zig-zag-coding

aString :: S.ByteString -> Builder
aString s = anInt (S.length s) <> byteString s

aCall :: Calls -> Builder
aCall (snps, indels) = mconcat
    [ aCommonCall snps, someGLs (vc_vars snps)
    , aCommonCall indels, someGLs (fst (vc_vars indels))
    , someIndels (snd (vc_vars indels)) ]
  where
    aCommonCall VarCall{..} = mconcat [ anInt vc_depth
                                      , anInt vc_mapq0
                                      , anInt vc_sum_mapq
                                      , anInt vc_sum_mapq2 ]

    someGLs gls = anInt (V.length gls) <> foldMap ((word8 . float2mini)

    someIndels is = anInt (length is) <> foldMap anIndel
    anIndel (d,is) = anInt d <> anInt (length is) <> foldMap (word8 . unN)

