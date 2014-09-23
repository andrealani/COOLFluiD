LIST ( APPEND MeshFEMMove_files
BasePrepare.cxx
BasePrepare.hh
FEMMove.cxx
FEMMove.hh
FEMMoveData.cxx
FEMMoveData.hh
ImposedMovementPrepare.cxx
ImposedMovementPrepare.hh
MeshFEMMove.hh
OscillatingAirfoilPrepare.cxx
OscillatingAirfoilPrepare.hh
StdPrepare.cxx
StdPrepare.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
UpdateMesh.cxx
UpdateMesh.hh
Valve3DCyclePrepare.cxx
Valve3DCyclePrepare.hh
ValveCyclePrepare.cxx
ValveCyclePrepare.hh
ValveCyclePrepare2.cxx
ValveCyclePrepare2.hh
ValveMovementPrepare.cxx
ValveMovementPrepare.hh
)

LIST ( APPEND MeshFEMMove_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( MeshFEMMove )
CF_WARN_ORPHAN_FILES()
