IF (CF_HAVE_PARDISO)

LIST ( APPEND Pardiso_files
  Pardiso.cxx
  PardisoData.cxx
  PardisoData.hh
  Pardiso.hh
  PardisoMatrix.cxx
  PardisoMatrix.hh
  PardisoModule.hh
  PardisoVector.cxx
  PardisoVector.hh
  StdSetup.cxx
  StdSetup.hh
  StdSolveSys.cxx
  StdSolveSys.hh
  StdUnSetup.cxx
  StdUnSetup.hh
)

LIST ( APPEND Pardiso_includedirs ${CF_PARDISO_INCLUDE} )
LIST ( APPEND Pardiso_libs ${CF_PARDISO_LIB})
LIST ( APPEND Pardiso_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( Pardiso )
CF_WARN_ORPHAN_FILES()
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/PardisoConfig.hh.in ${CMAKE_CURRENT_SOURCE_DIR}/PardisoConfig.hh)

ENDIF (CF_HAVE_PARDISO)


