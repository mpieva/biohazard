set -e
haddock=""

run_with() {
	echo -e "\e[91mTesting with GHC $1\e[0m"
	cabal install -w ghc-$1 --only-dep --force-reinstalls --disable-library-profiling --disable-profiling -O0
	cabal configure -w ghc-$1 --enable-tests --disable-library-profiling --disable-profiling -O0
	cabal build
	cabal test
	echo -e "\e[32mDone with GHC $1\e[0m"
	if [ -z "$haddock" ] ; then
		echo -e "\e[91mTesting Haddock\e[0m"
	       	cabal haddock
		haddock="done"
		echo -e "\e[32mDone with Haddock\e[0m"
	fi
	cabal clean
}

ghc_versions() {
        runghc << EOF
import Distribution.Compiler
import Distribution.PackageDescription
import Distribution.PackageDescription.Parse
import Distribution.Verbosity
import Distribution.Version

main = do
        gpd <- readGenericPackageDescription normal "biohazard.cabal"
        mapM_ (putStrLn . showVersion) $ concat
                [ foldVersionRange [] return (const []) (const []) (++) (\_ _ -> []) vr
                | (GHC, vr) <- testedWith (packageDescription gpd) ]
EOF
}

cabal sandbox init
for v in `ghc_versions` ; do
        run_with $v
done
