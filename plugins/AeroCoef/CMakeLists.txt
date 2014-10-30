LIST ( APPEND AeroCoef_files
AeroCoef.hh
AeroCoefStopCondition.hh
AeroCoefStopCondition.cxx
)

LIST ( APPEND AeroCoef_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( AeroCoef )

##################################################################

LIST ( APPEND AeroCoefFS_files
AeroCoefFS.hh
NavierStokes2DConsComputeAero.hh
NavierStokes2DConsComputeAero.cxx
NavierStokes3DConsComputeAero.hh
NavierStokes3DConsComputeAero.cxx
NavierStokes2DConsComputeAeroHO.hh
NavierStokes2DConsComputeAeroHO.cxx
NavierStokes3DConsComputeAeroHO.hh
NavierStokes3DConsComputeAeroHO.cxx
NavierStokes2DConsComputeAeroHOIsoP2.hh
NavierStokes2DConsComputeAeroHOIsoP2.cxx
NavierStokesConsComputeAero.hh
NavierStokesConsComputeAero.cxx
NavierStokesSkinFrictionHeatFluxFS.hh
NavierStokesSkinFrictionHeatFluxFS.ci
NavierStokesSkinFrictionHeatFluxFS.cxx
)

LIST ( APPEND AeroCoefFS_requires_mods FluctSplit )
LIST ( APPEND AeroCoefFS_libs FluctSplit NavierStokes Framework )
CF_ADD_PLUGIN_LIBRARY ( AeroCoefFS )

##################################################################

LIST ( APPEND AeroCoefFVM_files
AeroCoefFVM.hh
AeroForcesFVMCC.hh
AeroForcesFVMCC.cxx
NavierStokesSkinFrictionHeatFluxCC.hh
NavierStokesSkinFrictionHeatFluxCC.cxx
NavierStokesSkinFrictionHeatFluxCC3D.hh
NavierStokesSkinFrictionHeatFluxCC3D.ci
NavierStokesSkinFrictionHeatFluxCC3D.cxx
NavierStokesBLExtractionCC.hh
NavierStokesBLExtractionCC.cxx
Euler2DConsComputeAeroFVMCC.hh
Euler2DConsComputeAeroFVMCC.cxx
Extract2DSectionCC.hh
Extract2DSectionCC.cxx
)

LIST ( APPEND AeroCoefFVM_requires_mods FiniteVolume )
LIST ( APPEND AeroCoefFVM_libs FiniteVolume NavierStokes Framework )
CF_ADD_PLUGIN_LIBRARY ( AeroCoefFVM )

##################################################################

LIST ( APPEND AeroCoefFVMNEQ_files
AeroCoefFVMNEQ.hh
NavierStokesSkinFrictionHeatFluxCCNEQ.hh
NavierStokesSkinFrictionHeatFluxCCNEQ.cxx
NavierStokesSkinFrictionHeatFluxCCNEQ3D.cxx
)

LIST ( APPEND AeroCoefFVMNEQ_requires_mods FiniteVolume NEQ )
LIST ( APPEND AeroCoefFVMNEQ_libs FiniteVolume NEQ AeroCoefFVM )
CF_ADD_PLUGIN_LIBRARY ( AeroCoefFVMNEQ )

##################################################################

LIST ( APPEND DataProcessingHeat_files
DataProcessingHeat.hh
Heat2DPrimComputeFlux.hh
Heat2DPrimComputeFlux.cxx
)

LIST ( APPEND DataProcessingHeat_requires_mods Heat )

LIST ( APPEND DataProcessingHeat_libs Heat )
CF_ADD_PLUGIN_LIBRARY ( DataProcessingHeat )

##################################################################

LIST ( APPEND AeroCoefSpectralFD_files
AeroCoefSpectralFD.hh
EulerConsComputeAeroSpectralFD.cxx
EulerConsComputeAeroSpectralFD.hh
NavierStokesConsComputeAeroSpectralFD.cxx
NavierStokesConsComputeAeroSpectralFD.hh
NavierStokesConsComputeAeroSpectralFDCompact.cxx
NavierStokesConsComputeAeroSpectralFDCompact.hh
)

LIST ( APPEND AeroCoefSpectralFD_requires_mods SpectralFD NavierStokes )
LIST ( APPEND AeroCoefSpectralFD_cflibs SpectralFD NavierStokes Framework )
CF_ADD_PLUGIN_LIBRARY ( AeroCoefSpectralFD )

##################################################################

LIST ( APPEND AeroCoefDG_files
AeroCoefDG.hh
EulerComputeAeroDG.hh
EulerComputeAeroDG.cxx
NavierStokesComputeAeroDG.hh
NavierStokesComputeAeroDG.cxx
)

LIST ( APPEND AeroCoefDG_requires_mods DiscontGalerkin NavierStokes )
LIST ( APPEND AeroCoefDG_cflibs DiscontGalerkin NavierStokes Framework )
CF_ADD_PLUGIN_LIBRARY ( AeroCoefDG )

##################################################################

LIST ( APPEND AeroCoefSpectralFV_files
AeroCoefSpectralFV.hh
EulerConsComputeAeroSpectralFV.cxx
EulerConsComputeAeroSpectralFV.hh
)

LIST ( APPEND AeroCoefSpectralFV_requires_mods SpectralFV NavierStokes )
LIST ( APPEND AeroCoefSpectralFV_libs SpectralFV NavierStokes Framework )

CF_ADD_PLUGIN_LIBRARY ( AeroCoefSpectralFV )

##################################################################
CF_WARN_ORPHAN_FILES()
