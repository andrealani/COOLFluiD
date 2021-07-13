LIST ( APPEND PoissonNEQ_files
PoissonNEQ.hh
PoissonNEQ2DRhoivtTvToCons.cxx
PoissonNEQ2DRhoivtTvToCons.hh
PoissonNEQConvVarSet.ci
PoissonNEQConvVarSet.cxx
PoissonNEQConvVarSet.hh
PoissonNEQDiffVarSet.ci
PoissonNEQDiffVarSet.cxx
PoissonNEQDiffVarSet.hh
PoissonNEQPhysicalModel.ci
PoissonNEQPhysicalModel.cxx
PoissonNEQPhysicalModel.hh
PoissonNEQTerm.ci
PoissonNEQTerm.hh
)

LIST ( APPEND PoissonNEQ_requires_mods NavierStokes NEQ )
LIST ( APPEND PoissonNEQ_cflibs Framework NavierStokes NEQ )
CF_ADD_PLUGIN_LIBRARY ( PoissonNEQ )

CF_WARN_ORPHAN_FILES()
#ADD_SUBDIRECTORY ( testcases )
