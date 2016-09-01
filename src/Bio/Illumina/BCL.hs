module Bio.Illumina.BCL where

-- ^ Handling of Illumina BCL files.
-- We will support plain BCL, gzipped BCL and bgzf'ed BCL.  Plain BCL
-- starts with a cluster count (4 bytes, little-endian).  Base calls
-- follow with one byte per base:  bits [0..1] encode the base in the
-- order ACGT, bits 2..7 contain the quality score.
--
-- We will have to read from many files, so reading reasonably sized
-- blocks is imperative.  The typical bcl file on a MiSeq is 0.5MB, on a
-- HiSeq it's about 3MB.  We simply read them completely---this requires
-- 0.5-1GB of memory on a typical run, which shouldn't be a problem.
-- It's more if decompression is necessary, but still reasonable.

-- The BCLs come with a companion 'filter' file.  These start with three
-- header words:  zero, format version number, number of clusters.  The
-- remainder is one byte(!) per cluster, bit 0 is the filter flag.

