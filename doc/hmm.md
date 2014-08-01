Suppose we want to read/write/(de-)compress bgzf files in parallel, but
(mostly) keep the iteratee abstraction...  how to do that?

We definitely have to read ahead.  If we implement this on top of seek
excptions, we have to discard the read ahead on seek, and lose
pipelining.  Otoh, a bgzf decoder has to catch seek exceptions anyway,
it might as well take care of the cleanup.

If we simply read ahead and 'spark' decompression of the blocks, is that
already good enough?  Even if decompression is a (safe!) foreign call?
Hard to tell, as it would involve some unsafePerformIO...  However,
assuming the threaded runtime, multiple foreign calls (to libz) can run
in parallel.  So maybe we simply start multiple (de-) compressor
threads.  The async ackage might help here?


Interface to compression codecs:  A size limited channel, each item we
plug in contains an MVAR to store the result in.  After (de-)
compression, the result is put into the MVar.

Interface to BGZF reader:  The reader is a separate thread, it receives
request on a channel and produces MVars containing chunks on a channel.
Requests can be seeking to a position (reading starts immediately and
continues indefinitely), cancellation (reading stops), reading of an
interval (the interval is read), termination (readng stops, the file is
closed).

The enumerator runs in the parent thread and communicates through
channels.  Ideally, we want to queue requests for intervals, but we also
want to support seeking and indefinite reading.  Does that mean we need
two communications channels?

