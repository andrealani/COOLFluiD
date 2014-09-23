IF (CF_HAVE_PLAS)
  LIST (APPEND PLaSTracking_includedirs ${PLAS_INCLUDE_DIR})
  LIST (APPEND PLaSTracking_libs ${PLAS_LIBRARY})
  LIST (APPEND PLaSTracking_cflibs Framework )
  LIST (APPEND PLaSTracking_files
    plasinterface.cxx
    PLaSModule.hh
    PLaSTracking.hh
    PLaSTracking.cxx
    PLaSTrackingData.hh
    PLaSTrackingData.cxx
    StdProcess.hh
    StdProcess.cxx
    StdSetup.hh
    StdSetup.cxx
    StdUnSetup.hh
    StdUnSetup.cxx
    StgImplementationStd.hh
    StgImplementationStd.cxx
    StgImplementationBlockable.hh
    StgImplementationBlockable.cxx
    StgImplementation.hh
  )
  CF_ADD_PLUGIN_LIBRARY ( PLaSTracking )
  CF_WARN_ORPHAN_FILES ()
ENDIF (CF_HAVE_PLAS)

