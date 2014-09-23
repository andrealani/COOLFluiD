LIST ( APPEND MeshTools_files
MeshTools.hh
ComputeMeshQuality.cxx
ComputeMeshQuality.hh
ComputeShortestDistance.cxx
ComputeShortestDistance.hh
ComputeWallDistance.cxx
ComputeWallDistance.hh
ComputeWallDistanceNewton.cxx
ComputeWallDistanceNewton.hh
ConcreteQualityCalculator.cxx
ConcreteQualityCalculator.hh
QualityCalculator.cxx
QualityCalculator.hh
ReadWallDistance.cxx
ReadWallDistance.hh
ChangeMesh.hh
ChangeMesh.cxx
)

LIST ( APPEND MeshTools_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( MeshTools )

##################################################################

LIST ( APPEND MeshToolsFVM_files
MeshToolsFVM.hh
ComputeWallDistanceVectorCC.cxx
ComputeWallDistanceVectorCC.hh
ComputeWallDistanceVector2CC.cxx
ComputeWallDistanceVector2CC.hh
ComputeWallDistanceNewtonCC.cxx
ComputeWallDistanceNewtonCC.hh
ComputeWallDistanceFVMCC.cxx
ComputeWallDistanceFVMCC.hh
)

IF (CF_HAVE_MPI)

  LIST ( APPEND MeshToolsFVM_files
  ComputeWallDistanceVector2CCMPI.cxx
  ComputeWallDistanceVector2CCMPI.hh
  )

ENDIF()

LIST ( APPEND MeshToolsFVM_requires_mods FiniteVolume )
LIST ( APPEND MeshToolsFVM_cflibs MeshTools FiniteVolume )

CF_ADD_PLUGIN_LIBRARY ( MeshToolsFVM )

##################################################################

CF_WARN_ORPHAN_FILES()
