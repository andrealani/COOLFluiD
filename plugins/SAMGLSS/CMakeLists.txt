IF (CF_HAVE_SAMG)

  LIST ( APPEND SAMGLSS_files
    SAMGLSS.hh
    SAMGLSS.cxx
    SAMGLSSData.hh
    SAMGLSSData.cxx
    SAMGLSSMatrix.hh
    SAMGLSSMatrix.cxx
    SAMGLSSVector.hh
    SAMGLSSVector.cxx
    SAMGLSSModule.hh
    StdSetup.hh
    StdSetup.cxx
    StdSolveSys.hh
    StdSolveSys.cxx
    StdUnSetup.hh
    StdUnSetup.cxx
    UserDefinitions.hh
    UserDefinitions.cxx
  )

  LIST ( APPEND SAMGLSS_includedirs ${SAMG_INCLUDE_DIR} )
  LIST ( APPEND SAMGLSS_libs ${SAMG_LIBRARY} )
  LIST ( APPEND SAMGLSS_cflibs Framework )
  ADD_DEFINITIONS (${SAMG_OPTIONS})

  CF_ADD_PLUGIN_LIBRARY ( SAMGLSS )

  CF_WARN_ORPHAN_FILES()

ENDIF (CF_HAVE_SAMG)

