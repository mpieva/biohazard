module BCL where

-- ^ Handling of Illumina BCL files.
-- We will support plain BCL and gzipped BCL.  Plain BCL starts with a
-- cluster count (4 bytes, little-endian).  Base calls follow with one
-- byte per base:  bits [0..1] encode the base in the order ACGT, bits
-- 2..7 contain the quality score.
--
-- When creating any other format, we need to combine with positions
-- from a pos, locs or clocs file.  Pos files are text, the first two
-- word on each line are x and y coordinate (signed decimal floating
-- point numbers).  They are rounded to integers and clamped to a
-- minimum of zero.
--
-- Locs files have three header word (4 bytes, little endian), the third
-- is the number of clusters.  They are followed by two floats (IEEE
-- single precision, little endian) for (x,y) of each cluster.
--
-- Clocs are complicated---need to find the Illumina docu.
