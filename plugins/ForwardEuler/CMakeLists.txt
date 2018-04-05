IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "ForwardEuler" )
ELSE ()
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

LIST ( APPEND ForwardEuler_files
CopySol.cxx
CopySol.hh
ForwardEuler.hh
ForwardEulerAPI.hh
FwdEuler.cxx
FwdEuler.hh
FwdEulerData.cxx
FwdEulerData.hh
StdPrepare.cxx
StdPrepare.hh
FSHOPrepare.cxx
FSHOPrepare.hh
StdSetup.cxx
StdSetup.hh
FSHOSetup.cxx
FSHOSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
FSHOUnSetup.cxx
FSHOUnSetup.hh
TwoLayerPrepare.cxx
TwoLayerPrepare.hh
TwoLayerSetup.cxx
TwoLayerSetup.hh
TwoLayerUnSetup.cxx
TwoLayerUnSetup.hh
TwoLayerUpdateSol.cxx
TwoLayerUpdateSol.hh
UpdateSol.cxx
UpdateSol.hh
)

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND ForwardEuler_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( ForwardEuler )
ELSE()
FOREACH (AFILE ${ForwardEuler_files} )
LIST(APPEND coolfluid-solver_files ../../plugins/ForwardEuler/${AFILE} )
ENDFOREACH(AFILE)
CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()

CF_WARN_ORPHAN_FILES()
