import Control.Exception                    ( try, IOException )
import Distribution.PackageDescription      ( PackageDescription(..) )
import Distribution.Simple
import Distribution.Simple.InstallDirs      ( docdir, mandir, CopyDest (NoCopyDest) )
import Distribution.Simple.LocalBuildInfo   ( LocalBuildInfo(..), absoluteInstallDirs )
import Distribution.Simple.Program.Db       ( ProgramDb, lookupProgram )
import Distribution.Simple.Program.Run      ( runProgramInvocation, programInvocation, progInvokeCwd )
import Distribution.Simple.Program.Types    ( ConfiguredProgram, simpleProgram )
import Distribution.Simple.Setup            ( copyDest, copyVerbosity, fromFlag, installVerbosity, haddockVerbosity )
import Distribution.Simple.Utils
import Distribution.Verbosity               ( Verbosity, moreVerbose )
import System.Exit                          ( exitSuccess )
import System.FilePath                      ( splitDirectories, joinPath, takeExtension, replaceExtension, (</>) )

main :: IO ()
main = do
  defaultMainWithHooks $ simpleUserHooks
    { postCopy = \ _ flags pkg lbi ->
         installManpages pkg lbi (fromFlag $ copyVerbosity flags) (fromFlag $ copyDest flags)

    , postInst = \ _ flags pkg lbi ->
         installManpages pkg lbi (fromFlag $ installVerbosity flags) NoCopyDest

    , postHaddock = \ _ flags pkg lbi ->
         runPdflatex pkg lbi (fromFlag $ haddockVerbosity flags)

    , hookedPrograms = [ simpleProgram "pdflatex" ]
    }
  exitSuccess

installManpages :: PackageDescription -> LocalBuildInfo -> Verbosity -> CopyDest -> IO ()
installManpages pkg lbi verbosity copy = do
    installOrdinaryFiles verbosity (mandir (absoluteInstallDirs pkg lbi copy))
        [ ("man", joinPath mp) | ("man":mp) <- map splitDirectories $ extraSrcFiles pkg ]

    installOrdinaryFiles' verbosity (docdir (absoluteInstallDirs pkg lbi copy))
            [ (buildDir lbi </> "latex", replaceExtension (last p) "pdf")
            | ("doc":p@(_:_)) <- map splitDirectories $ extraSrcFiles pkg
            , takeExtension (last p) == ".tex" ]

installOrdinaryFiles' :: Verbosity -> FilePath -> [(FilePath, FilePath)] -> IO ()
installOrdinaryFiles' verb dest = mapM_ go
  where
    go :: (FilePath, FilePath) -> IO (Either IOException ())
    go (base,src) = try $ installOrdinaryFile verb (base </> src) (dest </> src)

withLatex :: LocalBuildInfo -> (ConfiguredProgram -> IO ()) -> IO ()
withLatex lbi k = maybe (return ()) k $ lookupProgram (simpleProgram "pdflatex") $ withPrograms lbi

runPdflatex :: PackageDescription -> LocalBuildInfo -> Verbosity -> IO ()
runPdflatex pkg lbi verb =
    withLatex lbi $ \cmd -> do
        createDirectoryIfMissingVerbose verb True (buildDir lbi </> "latex")
        sequence_ [ runProgramInvocation (moreVerbose verb) $
                        (programInvocation cmd [ "-interaction=nonstopmode", ddir </> joinPath ("doc":f) ])
                        { progInvokeCwd = Just bdir }
                  | ("doc":f@(_:_)) <- map splitDirectories $ extraSrcFiles pkg
                  , takeExtension (last f) == ".tex" ]
  where
    bdir = buildDir lbi </> "latex"
    ddir = joinPath (map (const "..") $ splitDirectories bdir)
