LIST ( APPEND EmptyConvergenceMethod_files
EmptyConvergenceMethod.hh
EmptyConvergenceMethodAPI.hh
EmptyIterator.cxx
EmptyIterator.hh
EmptyIteratorData.cxx
EmptyIteratorData.hh
StdSetup.cxx
StdSetup.hh
)

LIST ( APPEND EmptyConvergenceMethod_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( EmptyConvergenceMethod )
CF_WARN_ORPHAN_FILES()
