LIST ( APPEND Poisson_files
Poisson.hh
PoissonConvCons.cxx
PoissonConvCons.hh
PoissonConvVarSet.cxx
PoissonConvVarSet.hh
PoissonModel.ci
PoissonModel.cxx
PoissonModel.hh
PoissonConvTerm.cxx
PoissonConvTerm.hh
PoissonDiffCons.cxx
PoissonDiffCons.hh
PoissonDiffTerm.cxx
PoissonDiffTerm.hh
PoissonDiffVarSet.cxx
PoissonDiffVarSet.hh
PoissonReacTerm.cxx
PoissonReacTerm.hh
)

LIST ( APPEND Poisson_cflibs Common MathTools Framework )

CF_ADD_PLUGIN_LIBRARY ( Poisson )

CF_WARN_ORPHAN_FILES()
#ADD_SUBDIRECTORY ( testcases )
