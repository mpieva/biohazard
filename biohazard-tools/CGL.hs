module CGL where

-- "Compact Genotype Likelihoods"
--
-- Serialize the results from genotype calling in a sensible way.  A
-- file consists of 'Block's, each 'Block' comes with a header
-- identifying the reference sequence, and many 'Variant' records.
-- The file itself gets a minimal header to facilitate versioning, but
-- not much else.  We might want to compress 'Block's.


encodeCGL :: Monad m => String -> Refs -> Enumeratee [EitherCall] S.ByteString m a
encodeCGL = undefined
    -- encode a header (key and the string)
    -- loop
    --   encode refseq and start
    --   collect calls and keep encoding
    --   if the block is too big or we get EOF, flush
    --   maybe compress, needs testing

