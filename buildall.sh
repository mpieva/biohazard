set -e
version=$( awk '/Version:/ { print $2 }' biohazard.cabal )
echo "Building and installing biohazard-${version}"

cabal configure --ghc-options=-Wall -O2 --global
if which cabal-src-install ; then
    cabal-src-install --global
else
    cabal install --global
fi
cabal sdist

echo "SCP to bioinf..."
scp "dist/biohazard-${version}.tar.gz" bioinf:htdocs/biohazard/
