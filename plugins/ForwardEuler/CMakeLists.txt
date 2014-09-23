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

LIST ( APPEND ForwardEuler_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( ForwardEuler )
CF_WARN_ORPHAN_FILES()
