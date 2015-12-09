module Bio.Pipes.Bgzf where

-- ^ Bgzf decompression for Pipes.


-- | We need (semi-) random access, and we may need to handle
-- end-of-file, though most the time we don't actually care.  It seems,
-- a 'Parser' is overkill.  Instead, reading a file is a 'Server',
-- un-BGZF is a Proxy.  BGZF is a Pipe, and writing is a Sink.

-- | Random access to a handle.  Do we terminate at all, if so, what do
-- we return?  (We may need to use this in a 'Parser', but then there
-- isn't much point in having a 'Server'.)  If we go with this
-- interface, it may need to be duplicated for the case where seeking
-- isn't desired or can't be supported.
readHandle :: Server FileOffset ByteString m r

-- | Random access to a handle.  More explicit version.  If this is
-- incorporated into a 'Parser', termination can be detected, reading
-- can be restarted from somewhere else.  The handle needs to be closed
-- by some higher authority.
readHandle :: FileOffset -> Producer ByteString m r

-- | BGZF decompression.  Same nonsense applies, so maybe this is the
-- right thing to do.
inflateBGZF :: Parser ByteString (Producer Block m) r


