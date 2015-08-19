LIST ( APPEND Poisson_files
Poisson.hh
PoissonConv3DCons.cxx
PoissonConv3DCons.hh
PoissonConv3DVarSet.cxx
PoissonConv3DVarSet.hh
PoissonModel.ci
PoissonModel.cxx
PoissonModel.hh
PoissonConvTerm.cxx
PoissonConvTerm.hh
PoissonConvVarSet.hh
PoissonDiff3DCons.cxx
PoissonDiff3DCons.hh
PoissonDiffTerm.cxx
PoissonDiffTerm.hh
PoissonDiffVarSet.cxx
PoissonDiffVarSet.hh
PoissonReacTerm.cxx
PoissonReacTerm.hh
)

#LIST ( APPEND ArcJet_requires_mods NavierStokes LTE )
LIST ( APPEND Poisson_cflibs Common MathTools Framework )

CF_ADD_PLUGIN_LIBRARY ( Poisson )

CF_WARN_ORPHAN_FILES()
#ADD_SUBDIRECTORY ( testcases )
