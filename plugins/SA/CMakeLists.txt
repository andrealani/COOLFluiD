LIST ( APPEND SA_files
EulerSAConsVarSet.ci
EulerSAConsVarSet.cxx
EulerSAConsVarSet.hh
EulerSAPvtToCons.ci
EulerSAPvtToCons.cxx
EulerSAPvtToCons.hh
EulerSAVarSet.ci
EulerSAVarSet.cxx
EulerSAVarSet.hh
NavierStokesSAPvt.hh
NavierStokesSAPvt.ci
NavierStokesSAPvt.cxx
NavierStokesSAVarSet.ci
NavierStokesSAVarSet.hh
Euler2DSAConsToPuvtInPuvt.cxx
Euler2DSAConsToPuvtInPuvt.hh
Euler2DSAConsToRoe.cxx
Euler2DSAConsToRoe.hh
Euler2DSALinearCons.cxx
Euler2DSALinearCons.hh
Euler2DSALinearRoe.cxx
Euler2DSALinearRoe.hh
Euler2DSAPuvtToConsInPuvt.cxx
Euler2DSAPuvtToConsInPuvt.hh
NavierStokesSACons.ci
NavierStokesSACons.cxx
NavierStokesSACons.hh
NavierStokesSAPhysicalModel.ci
NavierStokesSAPhysicalModel.cxx
NavierStokesSAPhysicalModel.hh
NavierStokesSAVarSetTypes.hh
SA.hh
)

LIST ( APPEND SA_cflibs Framework NavierStokes)
LIST ( APPEND SA_requires_mods NavierStokes)
CF_ADD_PLUGIN_LIBRARY ( SA )

LIST ( APPEND SALTE_files
SALTE.hh
EulerSALTEVarSet.cxx
EulerSAPvtLTEToCons.cxx
NavierStokesSAPvtLTE.ci
NavierStokesSAPvtLTE.cxx
NavierStokesSAPvtLTE.hh
)

LIST ( APPEND SALTE_cflibs Framework NavierStokes LTE SA )
LIST ( APPEND SALTE_requires_mods NavierStokes LTE SA )
CF_ADD_PLUGIN_LIBRARY ( SALTE )

CF_WARN_ORPHAN_FILES()
ADD_SUBDIRECTORY ( testcases )

