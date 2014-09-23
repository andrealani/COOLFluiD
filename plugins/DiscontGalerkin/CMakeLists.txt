LIST ( APPEND DiscontGalerkin_files

DiscontGalerkin.hh
DGElemTypeData.hh
DiscontGalerkinSolver.hh
DiscontGalerkinSolver.cxx
DiscontGalerkinSolverData.hh
DiscontGalerkinSolverData.cxx
StdSolveCells.hh
StdSolveCells.cxx
StdSolveFaces.hh
StdSolveFaces.cxx
StdUnSetup.cxx
StdUnSetup.hh
StdSetup.hh
StdSetup.cxx
InitState.hh
InitState.cxx
DG_MeshDataBuilder.hh
DG_MeshDataBuilder.cxx
StdBaseSolve.hh
StdBaseSolve.cxx
ReflexWallBC.hh
ReflexWallBC.cxx
WallBC.hh
WallBC.cxx
InletBC.hh
InletBC.cxx
SuperInletBC.hh
SuperInletBC.cxx
SuperOutletBC.hh
SuperOutletBC.cxx
DG_UpgradeMeshDataBuilder.hh
DG_UpgradeMeshDataBuilder.cxx
ContourDGGaussLegendre3LagrangeTetra.cxx
ContourDGGaussLegendre3LagrangeTetra.hh
ContourDGGaussLegendre5LagrangeTetra.cxx
ContourDGGaussLegendre5LagrangeTetra.hh
StdStabilization.hh
StdStabilization.cxx
ViscousSolveCells.hh
ViscousSolveCells.cxx
ViscousSolveFaces.hh
ViscousSolveFaces.cxx
ViscousBaseSolve.hh
ViscousBaseSolve.cxx
ViscousWallBC.hh
ViscousWallBC.cxx
ViscousInletBC.hh
ViscousInletBC.cxx
StdComputeResidual.hh
StdComputeResidual.cxx
)

CF_ADD_PLUGIN_LIBRARY ( DiscontGalerkin )

CF_WARN_ORPHAN_FILES()
