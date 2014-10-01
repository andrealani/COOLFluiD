LIST ( APPEND ICP_files
ICP.hh
ICPInductionConvVarSet.hh
ICPInductionConvVarSet.ci
ICPInductionConvVarSet.cxx
ICPInductionDiffVarSet.hh
ICPInductionDiffVarSet.ci
ICPInductionDiffVarSet.cxx
ICPLTE2DPuvtToCons.hh
ICPLTE2DPuvtToCons.cxx
ICPLTE2DdPuvtToPuvt.hh
ICPLTE2DdPuvtToPuvt.cxx
ICPLTE2DPuvtToConsInPuvt.cxx
ICPLTE2DPuvtToConsInPuvt.hh
ICPLTEPhysicalModel.hh
ICPLTEPhysicalModel.ci
ICPLTEPhysicalModel.cxx
ICPReactionTerm.hh
ICPReactionTerm.ci
)

LIST ( APPEND ICPNEQ_files
ICPNEQ.hh
ICPNEQ2DPvtToPivtTv.hh
ICPNEQ2DPvtToPivtTv.cxx
ICPNEQ2DPuvtToRhoivtTv.hh
ICPNEQ2DPuvtToRhoivtTv.cxx
ICPNEQ2DPivtTvToCons.hh
ICPNEQ2DPivtTvToCons.cxx
ICPNEQ2DRhoivtToCons.hh
ICPNEQ2DRhoivtToCons.cxx
ICPNEQ2DRhoivtTvToCons.hh
ICPNEQ2DRhoivtTvToCons.cxx
ICPNEQInductionDiffVarSet.cxx
ICPNEQInductionConvVarSet.cxx
ICPNEQPhysicalModel.hh
ICPNEQPhysicalModel.ci
ICPNEQPhysicalModel.cxx
)

LIST ( APPEND ICP_requires_mods NavierStokes )
LIST ( APPEND ICP_cflibs Framework )
LIST ( APPEND ICP_libs NavierStokes ICP LTE )

LIST ( APPEND ICPNEQ_requires_mods NavierStokes ICP NEQ )
LIST ( APPEND ICPNEQ_libs NavierStokes ICP NEQ )

CF_ADD_PLUGIN_LIBRARY ( ICP )
CF_ADD_PLUGIN_LIBRARY ( ICPNEQ )

CF_WARN_ORPHAN_FILES()
#ADD_SUBDIRECTORY ( testcases )
