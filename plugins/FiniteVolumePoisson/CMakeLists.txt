LIST ( APPEND FiniteVolumePoisson_files
CurrentCondition.cxx
CurrentCondition.hh
FiniteVolumePoisson.hh
NeumannCondition.cxx
NeumannCondition.hh
PureDiffFlux.cxx
PureDiffFlux.hh
)

LIST ( APPEND FiniteVolumePoisson_requires_mods Poisson FiniteVolume )
LIST ( APPEND FiniteVolumePoisson_cflibs Framework Poisson FiniteVolume ) 

CF_ADD_PLUGIN_LIBRARY ( FiniteVolumePoisson )

CF_WARN_ORPHAN_FILES()
