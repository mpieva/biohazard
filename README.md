biohazard
=========

[![Build Status](https://travis-ci.org/udo-stenzel/biohazard.svg?branch=master)](https://travis-ci.org/udo-stenzel/biohazard)

This is a library for Bioinformatics, mostly centered on the BAM file
format for Next Generation Sequencing data.  It comes with a small set
of command line tools for operations on BAM files:

* `bam-fixpair`: brings mates from paired end runs together and fixes
  their flags

* `bam-meld`: melds multiple BAM files together, retaining only the best
  alignment for each read
  
* `bam-rmdup`: removes PCR duplicates from BAM files and computes a
  consensus sequence for each cluster of dups

* `bam-rewrap`: wraps alignments around the origin of a circular reference.
  (Use this to get sensible alignments to mitochondria.  `bam-rmdup` 
  includes similar functionality.)

Installation
------------

Biohazard uses Cabal, the standard installation mechanism for Haskell.  
To install, follow these steps:

* install a useable Haskell environment, either
 * install the Haskell Platform or have someone install it (see http://haskell.org/platform/ for details), or
 * install GHC (see http://haskell.org/ghc) and bootstrap Cabal (see
http://haskell.org/cabal), and
 * run `cabal update` (takes a while to download the current package list),
* run `cabal install <pkg-list>` in the biohazard directory (takes even longer).

That's it.  The `cabal update` step will probably tell you that there's a new version of
`cabal-install`, which you can install.  Feel free to do that, it's not
necessary, but doesn't hurt either.  

Cabal will download and install a number of dependencies,
then install Biohazard.  Before running `cabal install`, you can run 
`cabal install --dry-run -v` to check what it will do.  It should tell 
you that it's going to install a number of new libraries, but it 
shouldn't try any updates.

When done, on an unmodified Cabal setup, you will find the binaries in 
`${HOME}/cabal/bin`.  Cabal can install them in a different place, please 
refer to the Cabal documentation at http://www.haskell.org/cabal/ if 
you need that.  Sometimes, repeated installations and re-installations can result 
in a thoroughly unusable state of the Cabal package collection.  If you get error 
messages that just don't make sense anymore, please refer to 
http://www.vex.net/~trebla/haskell/sicp.xhtml; among other useful things, it 
tells you how to wipe a package database without causing more destruction.
