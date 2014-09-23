LIST ( APPEND EmptySpaceMethod_files
Empty.hh
EmptyMeshDataBuilder.cxx
EmptyMeshDataBuilder.hh
EmptySolver.cxx
EmptySolver.hh
EmptySolverData.cxx
EmptySolverData.hh
EmptyStrategy.cxx
EmptyStrategy.hh
StdSetup.cxx
StdSetup.hh
StdSolve.cxx
StdSolve.hh
StdUnSetup.cxx
StdUnSetup.hh
)

LIST ( APPEND EmptySpaceMethod_cflibs Framework)
CF_ADD_PLUGIN_LIBRARY ( EmptySpaceMethod )
CF_WARN_ORPHAN_FILES()
