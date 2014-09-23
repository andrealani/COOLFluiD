LIST ( APPEND LagrangianSolver_files
LagrangianSolverModule.hh
LagrangianSolver.hh
ParallelVector/ParallelVector.hh
ParallelVector/ParallelVector.cxx
ParticleTracking/ParticleTracking.hh
ParticleTracking/ParticleTracking.cxx
ParticleTracking/ParticleTracking2D.hh
ParticleTracking/ParticleTracking2D.cxx
ParticleTracking/ParticleTracking3D.hh
ParticleTracking/ParticleTracking3D.cxx
ParticleTracking/ParticleTrackingAxi.hh
ParticleTracking/ParticleTrackingAxi.cxx
SendBuffer/SendBuffer.hh
SendBuffer/SendBuffer.cxx
)

LIST ( APPEND LagrangianSolver_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( LagrangianSolver )

CF_WARN_ORPHAN_FILES()
