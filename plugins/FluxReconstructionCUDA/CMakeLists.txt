IF(CF_HAVE_CUDA)
LIST ( APPEND FluxReconstructionCUDA_files
     FluxReconstructionCUDA.hh
     StencilCUDASetup.cxx
     StencilCUDASetup.hh
     ConvRHSFluxReconstructionCUDA.hh
     ConvRHSFluxReconstructionCUDA.ci
     ConvRHSFluxReconstructionCUDA.cu
     ConvDiffRHSFluxReconstructionCUDA.hh
     ConvDiffRHSFluxReconstructionCUDA.ci
     ConvDiffRHSFluxReconstructionCUDA.cu
     NSBndGradientComputerCUDA.hh
     NSBndGradientComputerCUDA.cxx
     DiffBndCorrectionsRHSFluxReconstructionNSCUDA.hh
     DiffBndCorrectionsRHSFluxReconstructionNSCUDA.cxx
     PrepareCUDA.hh
     PrepareCUDA.cxx
)
    
# StencilCUDASetup.cxx or some other DUMMY file is 
# needed in order to properly link this module

LIST ( APPEND FluxReconstructionCUDA_requires_mods NavierStokes FluxReconstructionMethod FluxReconstructionNavierStokes)
LIST ( APPEND FluxReconstructionCUDA_cflibs NavierStokes FluxReconstructionMethod FluxReconstructionNavierStokes)
LIST ( APPEND FluxReconstructionCUDA_includedirs ${MPI_INCLUDE_DIR} ${CUDA_INCLUDE_DIR} )
LIST ( APPEND FluxReconstructionCUDA_libs ${CUDA_LIBRARIES} ) 

IF (CF_HAVE_PARALUTION)
LIST ( APPEND FluxReconstructionCUDA_includedirs ${PARALUTION_INCLUDE_DIR} )
LIST ( APPEND FluxReconstructionCUDA_libs ${PARALUTION_LIBRARIES} )
ENDIF()

CF_ADD_PLUGIN_LIBRARY ( FluxReconstructionCUDA )

ENDIF()
