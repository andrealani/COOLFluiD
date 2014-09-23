IF (CF_HAVE_TRILINOS)

LIST ( APPEND Trilinos_files
StdParSetup.cxx
StdParSetup.hh
StdParSolveSys.cxx
StdParSolveSys.hh
StdParUnSetup.cxx
StdParUnSetup.hh
Trilinos.hh
TrilinosIdxMapping.cxx
TrilinosIdxMapping.hh
TrilinosLSS.cxx
TrilinosLSS.hh
TrilinosLSSData.cxx
TrilinosLSSData.hh
TrilinosMatrix.cxx
TrilinosMatrix.hh
TrilinosOptions.cxx
TrilinosOptions.hh
TrilinosVector.cxx
TrilinosVector.hh
)

LIST ( APPEND Trilinos_includedirs ${TRILINOS_INCLUDE_DIR} )
#LIST ( APPEND Trilinos_libs ${TRILINOS_LIBRARIES} ${CF_BLASLAPACK_LIBRARIES} )
LIST ( APPEND Trilinos_libs ${TRILINOS_LIBRARIES} )
LIST ( APPEND Trilinos_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( Trilinos )
CF_WARN_ORPHAN_FILES()

ENDIF (CF_HAVE_TRILINOS)

