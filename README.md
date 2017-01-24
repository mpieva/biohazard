biohazard
=========

This is a library for Bioinformatics, mostly centered on the BAM file
format for Next Generation Sequencing data.


Installation
------------

Biohazard uses Cabal, the standard installation mechanism for Haskell.  
To install, follow these steps:

* install a useable Haskell environment, either
 * install the Haskell Platform or have someone install it (see http://haskell.org/platform/ for details), or
 * install GHC (see http://haskell.org/ghc) and bootstrap Cabal (see
http://haskell.org/cabal), and
 * run `cabal update` (takes a while to download the current package list),
* run `cabal install
  https://bitbucket.org/ustenzel/biohazard/get/0.6.10.tar.gz`
  (takes even longer).

Sometimes, repeated installations and re-installations can result in a
thoroughly unusable state of the Cabal package collection.  If you get
error messages that just don't make sense anymore, please refer to
http://www.vex.net/~trebla/haskell/sicp.xhtml; among other useful
things, it tells you how to wipe a package database without causing more
destruction.
