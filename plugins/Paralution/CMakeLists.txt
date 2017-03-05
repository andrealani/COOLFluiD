IF (CF_HAVE_PARALUTION)


LIST ( APPEND ParalutionI_files
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
Paralution.hh
ParalutionIdxMapping.cxx
ParalutionIdxMapping.hh
ParalutionLSS.cxx
ParalutionLSS.hh
ParalutionLSSData.cxx
ParalutionLSSData.hh
ParalutionMatrix.cxx
ParalutionMatrix.hh
ParalutionVector.cxx
ParalutionVector.hh
StdSolveSys.cxx
StdSolveSys.hh
)



LIST ( APPEND ParalutionI_includedirs ${MPI_INCLUDE_DIR} ${PARALUTION_INCLUDE_DIR} )

LIST ( APPEND ParalutionI_cflibs Common Framework ${PARALUTION_LIBRARIES} )


IF(CF_HAVE_CUDA)
LIST ( APPEND ParalutionI_files
    ParalutionMatrix.cu
)
LIST ( APPEND ParalutionI_includedirs ${CUDA_INCLUDE_DIR} )
LIST ( APPEND ParalutionI_cflibs ${CUDA_LIBRARIES} )
ENDIF(CF_HAVE_CUDA)

CF_ADD_PLUGIN_LIBRARY ( ParalutionI )

CF_WARN_ORPHAN_FILES()

ENDIF (CF_HAVE_PARALUTION)

