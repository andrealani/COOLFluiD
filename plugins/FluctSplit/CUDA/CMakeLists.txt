IF(CF_HAVE_CUDA)
LIST ( APPEND FluctSplitCUDA_files
     FluctSplitCUDA.hh
     LDAC_CUDA.hh
     LDAC_CUDA.cu
     LDASchemeCSysCUDA.hh
     LDASchemeCSysCUDA.cxx
#     ComputeRHSCUDA.hh
#     ComputeRHSCUDA.cxx
)

LIST ( APPEND FluctSplitCUDA_requires_mods FluctSplit )
LIST ( APPEND FluctSplitCUDA_cflibs FluctSplit )
LIST ( APPEND FluctSplitCUDA_includedirs ${CUDA_INCLUDE_DIR} )
LIST ( APPEND FluctSplitCUDA_libs ${CUDA_LIBRARIES} )

CF_ADD_PLUGIN_LIBRARY ( FluctSplitCUDA )
ENDIF()
