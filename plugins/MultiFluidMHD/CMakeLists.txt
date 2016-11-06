LIST ( APPEND MultiFluidMHD_files
MultiFluidMHD.hh
EulerMFMHDTerm.cxx
EulerMFMHDTerm.hh
DiffMFMHD2DHalfRhoiViTi.cxx
DiffMFMHD2DHalfRhoiViTi.hh
DiffMFMHD2DHalfVarSet.cxx
DiffMFMHD2DHalfVarSet.hh
DiffMFMHD2DRhoiViTi.hh
DiffMFMHD2DRhoiViTi.cxx
DiffMFMHD2DVarSet.hh
DiffMFMHD2DVarSet.cxx
DiffMFMHD3DRhoiViTi.cxx
DiffMFMHD3DRhoiViTi.hh
DiffMFMHD3DVarSet.hh
DiffMFMHD3DVarSet.cxx
DiffMFMHDTerm.cxx
DiffMFMHDTerm.hh
DiffMFMHDVarSet.cxx
DiffMFMHDVarSet.hh
Euler2DHalfMFMHDConsToRhoiViTiInRhoiViTi.cxx
Euler2DHalfMFMHDConsToRhoiViTiInRhoiViTi.hh
Euler2DMFMHDConsToRhoiViTiInRhoiViTi.cxx
Euler2DMFMHDConsToRhoiViTiInRhoiViTi.hh
EulerMFMHD2DCons.cxx
EulerMFMHD2DCons.hh
EulerMFMHD2DRhoiViTi.cxx
EulerMFMHD2DRhoiViTi.hh
EulerMFMHD2DConsToRhoiViTi.cxx
EulerMFMHD2DConsToRhoiViTi.hh
EulerMFMHD2DHalfCons.cxx
EulerMFMHD2DHalfCons.hh
EulerMFMHD2DHalfConsToRhoiViTi.cxx
EulerMFMHD2DHalfConsToRhoiViTi.hh
EulerMFMHD2DHalfRhoiViTi.cxx
EulerMFMHD2DHalfRhoiViTi.hh
EulerMFMHD2DHalfRhoiViTiToCons.cxx
EulerMFMHD2DHalfRhoiViTiToCons.hh
EulerMFMHD2DRhoiViTiToCons.cxx
EulerMFMHD2DRhoiViTiToCons.hh
EulerMFMHD3DRhoiViTi.cxx
EulerMFMHD3DRhoiViTi.hh
EulerMFMHD3DRhoiViTiToCons.cxx
EulerMFMHD3DRhoiViTiToCons.hh
EulerMFMHD3DConsToRhoiViTi.cxx
EulerMFMHD3DConsToRhoiViTi.hh
Euler3DMFMHDConsToRhoiViTiInRhoiViTi.cxx
Euler3DMFMHDConsToRhoiViTiInRhoiViTi.hh
EulerMFMHD3DCons.cxx
EulerMFMHD3DCons.hh
MultiFluidMHDModel.ci
MultiFluidMHDModel.cxx
MultiFluidMHDModel.hh
MultiFluidMHDVarSet.ci
MultiFluidMHDVarSet.hh
ReacMFMHDTerm.cxx
ReacMFMHDTerm.hh
)


IF (CF_HAVE_CUDA) 
   LIST ( APPEND MultiFluidMHD_files
MultiFluidMHD2DHalfTwoSpeciesVarSetT.hh 
EulerMFMHD2DHalfConsT.hh                                                                                        
EulerMFMHD2DHalfConsToRhoiViTiT.hh                                                                            
EulerMFMHD2DHalfRhoiViTiT.hh                                                                                    
EulerMFMHD2DHalfRhoiViTiToConsT.hh  
#MultiFluidMHD2DTwoSpeciesVarSetT.hh 
EulerMFMHD2DConsT.hh                                                                                        
EulerMFMHD2DConsToRhoiViTiT.hh                                                                            
EulerMFMHD2DRhoiViTiT.hh                                                                                    
EulerMFMHD2DRhoiViTiToConsT.hh  
) 
ENDIF ()



LIST ( APPEND MultiFluidMHD_cflibs Framework Maxwell )
LIST ( APPEND MultiFluidMHD_requires_mods Maxwell )

CF_ADD_PLUGIN_LIBRARY ( MultiFluidMHD )
CF_WARN_ORPHAN_FILES()
ADD_SUBDIRECTORY ( testcases )
