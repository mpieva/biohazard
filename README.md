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

Installation
------------

Biohazard uses Cabal, the standard installation mechanism for Haskell.  
To install, follow these steps:

* install the Haskell Platform or have some one install it (see http://www.haskell.org/platform/ for details),
* run "cabal update" (takes a while to download the current package list),
* run "cabal install" in the biohazard directory (takes even longer).

That's it.  The "cabal update" step will tell you that there's a new version of
cabal-install, which you can install.  Feel free to do that, it's not
neccessary, but doesn't hurt either.  

"cabal install" will download and install a number of dependencies,
then install biohazard.  Before running "cabal install", you can run 
"cabal install --dry-run -v" to check what it will do.  It should tell 
you that it's going to install a number of new libraries, but it 
shouldn't try any updates.

When done, on an unmodified Cabal setup, you will find the binaries in 
${HOME}/cabal/bin.  Cabal can install them in a different place, please 
refer to the Cabal documentation at http://www.haskell.org/cabal/ if 
you need that.
