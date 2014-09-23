LIST ( APPEND MeshRigidMove_files
MeshRigidMove.hh
OscillatingAirfoil.cxx
OscillatingAirfoil.hh
Piston.cxx
Piston.hh
PistonConstantSpeed.cxx
PistonConstantSpeed.hh
PistonDynamics.cxx
PistonDynamics.hh
RigidMove.cxx
RigidMove.hh
RigidMoveData.cxx
RigidMoveData.hh
StdPrepare.cxx
StdPrepare.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
TestOscillation.cxx
TestOscillation.hh
UpdateMesh.cxx
UpdateMesh.hh
Wedge.cxx
Wedge.hh
)

LIST ( APPEND MeshRigidMove_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( MeshRigidMove )
CF_WARN_ORPHAN_FILES()
