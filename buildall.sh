set -e
version=$( awk '/Version:/ { print $2 }' biohazard.cabal )
echo "Building and installing biohazard-${version}"

cabal clean
cabal configure --ghc-options=-Wall -O2 --prefix=/home/public/usr64/stow/biohazard-${version}
cabal build
if which cabal-src-install ; then
    cabal-src-install --prefix=/home/public/usr64/stow/biohazard-${version}
else
    cabal install --prefix=/home/public/usr64/stow/biohazard-${version}
fi
cabal sdist

(
  cd /home/public/usr64/stow
  stow -D biohazard-[0-9]*
  stow biohazard-${version}
  tar -czf /var/tmp/biohazard-${version}_x86_64.tar.gz biohazard-${version}
)

echo "SCP to bioinf..."
scp dist/biohazard-${version}.tar.gz \
scp /var/tmp/biohazard-${version}_x86_64.tar.gz \
    bioinf:htdocs/biohazard/
rm /var/tmp/biohazard-${version}_x86_64.tar.gz 
