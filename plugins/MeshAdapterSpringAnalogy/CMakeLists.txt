LIST ( APPEND MeshAdapterSpringAnalogy_files
BallVertexCalculator.cxx
BallVertexCalculator.hh
BasePrepare.cxx
BasePrepare.hh
CoupledPrepare.cxx
CoupledPrepare.hh
ImposedMovementPrepare.cxx
ImposedMovementPrepare.hh
ImposedValvePrepare.cxx
ImposedValvePrepare.hh
MeshAdapterSpringAnalogy.hh
SpringAnalogy.cxx
SpringAnalogy.hh
SpringAnalogyData.cxx
SpringAnalogyData.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
UpdateMesh.cxx
UpdateMesh.hh
)

LIST ( APPEND MeshAdapterSpringAnalogy_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( MeshAdapterSpringAnalogy )
CF_WARN_ORPHAN_FILES()
