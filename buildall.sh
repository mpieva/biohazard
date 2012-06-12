set -e
version=$( awk '/Version:/ { print $2 }' biohazard.cabal )
echo "Building and installing biohazard-${version}"

cabal clean
cabal configure --ghc-options=-Wall -O2
cabal build
cabal install --force-reinstall
cabal sdist

echo "SCP to bioinf..."
scp dist/biohazard-${version}.tar.gz \
    bioinf:htdocs/biohazard/
