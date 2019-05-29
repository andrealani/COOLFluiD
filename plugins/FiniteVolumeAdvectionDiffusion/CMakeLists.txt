LIST ( APPEND FiniteVolumeAdvectionDiffusion_files
FiniteVolumeAdvectionDiffusion.hh
LinearAdvDiffFlux.cxx
)

LIST ( APPEND FiniteVolumeAdvectionDiffusion_requires_mods LinearAdv FiniteVolume )
LIST ( APPEND FiniteVolumeAdvectionDiffusion_cflibs Framework ShapeFunctions LinearAdv FiniteVolume ) 

CF_ADD_PLUGIN_LIBRARY ( FiniteVolumeAdvectionDiffusion )

CF_WARN_ORPHAN_FILES()
