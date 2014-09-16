import Distribution.Simple
import Distribution.Simple.Setup ( copyDest, copyVerbosity, fromFlag, installVerbosity )
import Distribution.PackageDescription      ( PackageDescription(..) )
import Distribution.Simple.LocalBuildInfo   ( LocalBuildInfo(..), absoluteInstallDirs )
import Distribution.Verbosity               ( Verbosity )
import Distribution.Simple.InstallDirs      ( mandir, CopyDest (NoCopyDest) )
import Distribution.Simple.Utils            ( installOrdinaryFiles  )
import System.FilePath                      ( splitDirectories, joinPath )
import System.Exit

main :: IO ()
main = do
  defaultMainWithHooks $ simpleUserHooks
    { postCopy = \ _ flags pkg lbi ->
         installManpages pkg lbi (fromFlag $ copyVerbosity flags)
              (fromFlag $ copyDest flags)
    , postInst = \ _ flags pkg lbi ->
         installManpages pkg lbi (fromFlag $ installVerbosity flags)
              NoCopyDest
    }
  exitSuccess

installManpages :: PackageDescription -> LocalBuildInfo -> Verbosity -> CopyDest -> IO ()
installManpages pkg lbi verbosity copy =
    installOrdinaryFiles verbosity (mandir (absoluteInstallDirs pkg lbi copy))
        [ ("man", joinPath mp) | ("man":mp) <- map splitDirectories $ extraSrcFiles pkg ]
