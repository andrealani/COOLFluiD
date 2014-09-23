LIST ( APPEND MeshLaplacianSmoothing_files
LaplacianSmoothing.cxx
LaplacianSmoothing.hh
LaplacianSmoothingData.cxx
LaplacianSmoothingData.hh
MeshLaplacianSmoothing.hh
StdPrepare.cxx
StdPrepare.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
UpdateMesh.cxx
UpdateMesh.hh
)

LIST ( APPEND MeshLaplacianSmoothing_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( MeshLaplacianSmoothing )
CF_WARN_ORPHAN_FILES()
