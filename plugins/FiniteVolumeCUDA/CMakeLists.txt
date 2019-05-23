IF(CF_HAVE_CUDA)
LIST ( APPEND FiniteVolumeCUDA_files
     FiniteVolumeCUDA.hh
     FVMCC_ComputeRHSCell.ci
     FVMCC_ComputeRHSCell.cu
     FVMCC_ComputeRHSCell.hh     
     FVMCC_ComputeRhsJacobCell.ci
     FVMCC_ComputeRhsJacobCell.cu
     FVMCC_ComputeRhsJacobCell.hh
     FVMCC_ComputeSourceRHSCell.ci
     FVMCC_ComputeSourceRHSCell.hh
     FVMCC_ComputeSourceRHSCellExec.cu
     FVMCC_ComputeSourceRHSCellExec.hh
     FVMCC_ComputeSourceRhsJacobCell.ci
     FVMCC_ComputeSourceRhsJacobCell.hh
     FVMCC_ComputeSourceRhsJacobCellExec.cu
     FVMCC_ComputeSourceRhsJacobCellExec.hh
     FVMCC_ComputeRHSCellMHD.cu
     StencilCUDASetup.cxx	
     StencilCUDASetup.hh
)
    
# StencilCUDASetup.cxx or some other DUMMY file is 
# needed in order to properly link this module

LIST ( APPEND FiniteVolumeCUDA_requires_mods MHD FiniteVolume FiniteVolumeMHD FiniteVolumeMaxwell Maxwell FiniteVolumeMultiFluidMHD MultiFluidMHD)
LIST ( APPEND FiniteVolumeCUDA_cflibs MHD FiniteVolume FiniteVolumeMHD FiniteVolumeMaxwell Maxwell FiniteVolumeMultiFluidMHD MultiFluidMHD)
#LIST ( APPEND FiniteVolumeCUDA_includedirs ${MPI_INCLUDE_DIR} ${CUDA_INCLUDE_DIR} )
#LIST ( APPEND FiniteVolumeCUDA_libs ${MPI_LIBRARIES} ${CUDA_LIBRARIES} ) 

CF_ADD_PLUGIN_LIBRARY ( FiniteVolumeCUDA )
ENDIF()

IF (CF_HAVE_CUDA AND CF_HAVE_PARALUTION)
LIST ( APPEND FiniteVolumeCUDAParalution_files
     DummyCom.cxx
     DummyCom.hh
     FiniteVolumeCUDAParalution.hh
     FVMCC_ComputeSourceRhsJacobCellParalution.cu
     FVMCC_ComputeSourceRhsJacobCellParalution.hh
)

LIST ( APPEND FiniteVolumeCUDAParalution_requires_mods ParalutionI MHD FiniteVolume FiniteVolumeMHD FiniteVolumeMaxwell Maxwell FiniteVolumeMultiFluidMHD MultiFluidMHD FiniteVolumeCUDA )
LIST ( APPEND FiniteVolumeCUDAParalution_cflibs ParalutionI MHD FiniteVolume FiniteVolumeMHD FiniteVolumeMaxwell Maxwell FiniteVolumeMultiFluidMHD MultiFluidMHD FiniteVolumeCUDA )
LIST ( APPEND FiniteVolumeCUDAParalution_includedirs ${PARALUTION_INCLUDE_DIR} )
LIST ( APPEND FiniteVolumeCUDAParalution_libs ${PARALUTION_LIBRARIES} )

CF_ADD_PLUGIN_LIBRARY ( FiniteVolumeCUDAParalution )
ENDIF()

