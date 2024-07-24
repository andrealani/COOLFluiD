LIST ( APPEND HyperPoisson_files
HyperPoisson.hh
HyperPoisson3DCons.cxx
HyperPoisson3DCons.hh
HyperPoisson3DVarSet.cxx
HyperPoisson3DVarSet.hh
HyperPoissonVarSet.hh
HyperPTerm.cxx
HyperPTerm.hh
HyperPoissonPhysicalModel.cxx
HyperPoissonPhysicalModel.hh
HyperPoissonPhysicalModel.ci
)

LIST ( APPEND HyperPoisson_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( HyperPoisson )
CF_WARN_ORPHAN_FILES()

#ADD_SUBDIRECTORY ( testcases )