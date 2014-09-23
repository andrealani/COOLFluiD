LIST ( APPEND LinEuler_files
CreateDampingCoefficient.hh
CreateDampingCoefficient.cxx
CreateMeanFlowAnalytic.hh
CreateMeanFlowAnalytic.cxx
LinearizedEuler.hh
LinEulerPhysicalModel.ci
LinEulerPhysicalModel.cxx
LinEulerPhysicalModel.hh
LinEulerTerm.cxx
LinEulerTerm.hh
LinEulerVarSet.hh
LinEuler2DCons.cxx
LinEuler2DCons.hh
LinEuler2DLinearCons.cxx
LinEuler2DLinearCons.hh
LinEuler2DVarSet.cxx
LinEuler2DVarSet.hh
LinEuler3DCons.cxx
LinEuler3DCons.hh
LinEuler3DLinearCons.cxx
LinEuler3DLinearCons.hh
LinEuler3DVarSet.cxx
LinEuler3DVarSet.hh
ReadMeanFlowIncomp.cxx
ReadMeanFlowIncomp.hh
ReadMeanFlowComp.cxx
ReadMeanFlowComp.hh
)

LIST ( APPEND LinEuler_libs Framework )

CF_ADD_PLUGIN_LIBRARY ( LinEuler )
CF_WARN_ORPHAN_FILES()
ADD_SUBDIRECTORY ( testcases )
