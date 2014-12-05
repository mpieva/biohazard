import Distribution.PackageDescription      ( PackageDescription(..) )
import Distribution.Simple
import Distribution.Simple.InstallDirs      ( docdir, mandir, CopyDest (NoCopyDest) )
import Distribution.Simple.LocalBuildInfo   ( LocalBuildInfo(..), absoluteInstallDirs )
import Distribution.Simple.Program          ( runProgram )
import Distribution.Simple.Program.Db       ( ProgramDb, lookupProgram )
import Distribution.Simple.Program.Types    ( ConfiguredProgram, simpleProgram )
import Distribution.Simple.Setup            ( copyDest, copyVerbosity, fromFlag, installVerbosity, buildVerbosity )
import Distribution.Simple.Utils            ( installOrdinaryFiles  )
import Distribution.Verbosity               ( Verbosity, moreVerbose )
import System.Exit                          ( exitSuccess )
import System.FilePath                      ( splitDirectories, joinPath, takeExtension, replaceExtension, (</>) )
import System.Directory                     ( getCurrentDirectory, setCurrentDirectory, createDirectoryIfMissing )

main :: IO ()
main = do
  defaultMainWithHooks $ simpleUserHooks
    { postCopy = \ _ flags pkg lbi ->
         installManpages pkg lbi (fromFlag $ copyVerbosity flags) (fromFlag $ copyDest flags)

    , postInst = \ _ flags pkg lbi ->
         installManpages pkg lbi (fromFlag $ installVerbosity flags) NoCopyDest

    -- , postBuild = \ _ flags pkg lbi ->
         -- runPdflatex pkg lbi (fromFlag $ buildVerbosity flags)

    , hookedPrograms = [ simpleProgram "pdflatex" ]
    }
  exitSuccess

installManpages :: PackageDescription -> LocalBuildInfo -> Verbosity -> CopyDest -> IO ()
installManpages pkg lbi verbosity copy = do
    installOrdinaryFiles verbosity (mandir (absoluteInstallDirs pkg lbi copy))
        [ ("man", joinPath mp) | ("man":mp) <- map splitDirectories $ extraSrcFiles pkg ]

    -- installOrdinaryFiles verbosity (docdir (absoluteInstallDirs pkg lbi copy))
            -- [ (buildDir lbi </> "latex", replaceExtension (last p) "pdf")
            -- | ("doc":p@(_:_)) <- map splitDirectories $ extraSrcFiles pkg
            -- , takeExtension (last p) == ".tex" ]

withLatex :: LocalBuildInfo -> (ConfiguredProgram -> IO ()) -> IO ()
withLatex lbi k = maybe (return ()) k $ lookupProgram (simpleProgram "pdflatex") $ withPrograms lbi

runPdflatex :: PackageDescription -> LocalBuildInfo -> Verbosity -> IO ()
runPdflatex pkg lbi verb =
    withLatex lbi $ \cmd -> do
        cwd <- getCurrentDirectory
        createDirectoryIfMissing True (buildDir lbi </> "latex")
        setCurrentDirectory (buildDir lbi </> "latex")
        sequence_ [ runProgram (moreVerbose verb) cmd [ "-interaction=batchmode", cwd </> joinPath ("doc":f) ]
                  | ("doc":f@(_:_)) <- map splitDirectories $ extraSrcFiles pkg
                  , takeExtension (last f) == ".tex" ]
        setCurrentDirectory cwd

