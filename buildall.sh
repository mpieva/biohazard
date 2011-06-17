set -e
version=$( awk '/Version:/ { print $2 }' biohazard.cabal )
echo "Building and installing biohazard-${version}"

cabal clean
cabal configure --ghc-options=-Wall
cabal build
cabal install
cabal sdist

echo "SCP to bioinf..."
scp dist/biohazard-${version}.tar.gz \
    bioinf:htdocs/biohazard/
