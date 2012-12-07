import Distribution.Simple
import Distribution.Simple.Setup ( copyDest, copyVerbosity, fromFlag, installVerbosity )
import Distribution.PackageDescription      ( PackageDescription(..) )
import Distribution.Simple.LocalBuildInfo   ( LocalBuildInfo(..), absoluteInstallDirs )
import Distribution.Verbosity               ( Verbosity ) 
import Distribution.Simple.InstallDirs      ( mandir, CopyDest (NoCopyDest) )
import Distribution.Simple.Utils            ( installOrdinaryFiles  )
import System.FilePath                      ( (</>) )
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
  exitWith ExitSuccess

manpages :: [FilePath]
manpages = map ("man1" </>) ["bam-rmdup.1"]

manDir :: FilePath
manDir = "man"

installManpages :: PackageDescription -> LocalBuildInfo -> Verbosity -> CopyDest -> IO ()
installManpages pkg lbi verbosity copy =
    installOrdinaryFiles verbosity
            (mandir (absoluteInstallDirs pkg lbi copy))
            (zip (repeat manDir) manpages)
