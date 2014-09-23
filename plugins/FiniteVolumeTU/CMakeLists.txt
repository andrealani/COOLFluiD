LIST ( APPEND FiniteVolumeTU_files
AUSMTU.cxx
FiniteVolumeTU.hh
Rotation3DSourceTerm.hh
Rotation3DSourceTerm.cxx
SubOutletEulerPvtTU.cxx
)

LIST ( APPEND FiniteVolumeTU_cflibs FiniteVolume FiniteVolumeNavierStokes NavierStokes )
LIST ( APPEND FiniteVolumeTU_requires_mods FiniteVolume FiniteVolumeNavierStokes NavierStokes )

CF_ADD_PLUGIN_LIBRARY ( FiniteVolumeTU )
CF_WARN_ORPHAN_FILES()
