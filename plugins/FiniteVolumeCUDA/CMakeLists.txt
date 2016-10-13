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
     FVMCC_ComputeSourceRHSCell.cu
     FVMCC_ComputeSourceRHSCell.hh     
     FVMCC_ComputeSourceRhsJacobCell.ci
     FVMCC_ComputeSourceRhsJacobCell.cu
     FVMCC_ComputeSourceRhsJacobCell.hh
     FVMCC_ComputeRHSCellMHD.cu
     StencilCUDASetup.cxx	
     StencilCUDASetup.hh
)
    
# StencilCUDASetup.cxx or some other DUMMY file is 
# needed in order to properly link this module

LIST ( APPEND FiniteVolumeCUDA_requires_mods MHD FiniteVolume FiniteVolumeMHD FiniteVolumeMaxwell Maxwell FiniteVolumeMultiFluidMHD MultiFluidMHD)
LIST ( APPEND FiniteVolumeCUDA_cflibs MHD FiniteVolume FiniteVolumeMHD FiniteVolumeMaxwell Maxwell FiniteVolumeMultiFluidMHD MultiFluidMHD)
LIST ( APPEND FiniteVolumeCUDA_includedirs ${MPI_INCLUDE_DIR} ${CUDA_INCLUDE_DIR} )
LIST ( APPEND FiniteVolumeCUDA_libs ${CUDA_LIBRARIES} )

CF_ADD_PLUGIN_LIBRARY ( FiniteVolumeCUDA )

ENDIF()
