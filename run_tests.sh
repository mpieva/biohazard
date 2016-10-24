set -e

cabal sandbox init

echo "Testing with GHC 7.8"
cabal configure -w ghc-7.8.4 --disable-library-profiling --disable-profiling -O0
cabal build
cabal test
echo "Done with GHC 7.8"

echo "Testing with GHC 7.10"
cabal configure -w ghc-7.10.1 --disable-library-profiling --disable-profiling -O0
cabal build
cabal test
echo "Done with GHC 7.10"

echo "Testing with GHC 8.0"
cabal configure -w ghc-8.0.1 --disable-library-profiling --disable-profiling -O0
cabal build
cabal test
cabal haddock
echo "Done with GHC 8.0"
