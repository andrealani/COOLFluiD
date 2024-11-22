LIST ( APPEND PoissonNEQ_files
PoissonNEQ.hh
PoissonNEQ2DRhoivtTvToCons.cxx #Vatsalya New addition
PoissonNEQ2DRhoivtTvToCons.hh   #Vatsalya New addition
PoissonNEQ3DRhoivtTvToCons.cxx  #Vatsalya New addition
PoissonNEQ3DRhoivtTvToCons.hh   #Vatsalya New addition
PoissonNEQ3DRhovtXiToRhoivtTv.cxx #Vatsalya: New addition
PoissonNEQ3DRhovtXiToRhoivtTv.hh # Vatsalya: New addition
PoissonNEQ2DRhoivtToCons.cxx
PoissonNEQ2DRhoivtToCons.hh
PoissonNEQ3DRhoivtToCons.cxx
PoissonNEQ3DRhoivtToCons.hh
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
