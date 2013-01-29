set -e
version=$( awk '/Version:/ { print $2 }' biohazard.cabal )
echo "Building and installing biohazard-${version}"

# STOW="/home/public/usr64/stow"
# DEST="${STOW}/biohazard-${version}"
STAR="dist/biohazard-${version}.tar.gz"
BTAR="dist/biohazard-${version}_x86_64.tar.gz"

cabal clean
cabal configure --ghc-options=-Wall -O2 --global
cabal build
if which cabal-src-install ; then
    cabal-src-install --global
else
    cabal install --global
fi
cabal sdist

# (
  # cd "${STOW}"
  # stow -D biohazard-[0-9]*
  # stow "biohazard-${version}"
# )
# tar -czf "${BTAR}" -C "${STOW}" "biohazard-${version}"

echo "SCP to bioinf..."
scp "${STAR}" bioinf:htdocs/biohazard/
