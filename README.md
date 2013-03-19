biohazard
=========

This is a library for Bioinformatics, mostly centered on the BAM file
format for Next Generation Sequencing Data.  It comes with a small set
of command line tools for operations on BAM files:

* bam-fixpair: brings mates from paired end runs together and fixes
  their flags

* bam-meld: melds multiple BAM files together, retaining only the best
  alignment for each read
  
* bam-rmdup: removes PCR duplicates from BAM files and computes a
  consensus sequence for each cluster of dups

