LIST ( APPEND FiniteVolumeRadiation_files
FiniteVolumeRadiation.hh
Radiation.cxx
Radiation.hh
)

LIST ( APPEND FiniteVolumeRadiation_cflibs NavierStokes LTE FiniteVolume FiniteVolumeNavierStokes )
LIST ( APPEND FiniteVolumeRadiation_requires_mods NavierStokes LTE FiniteVolume FiniteVolumeNavierStokes )

CF_ADD_PLUGIN_LIBRARY ( FiniteVolumeRadiation )

CF_WARN_ORPHAN_FILES()
