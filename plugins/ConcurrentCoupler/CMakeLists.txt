LIST ( APPEND ConcurrentCoupler_files
ConcurrentCoupler.hh
ConcurrentCouplerMethod.cxx
ConcurrentCouplerMethod.hh
ConcurrentCouplerData.cxx
ConcurrentCouplerData.hh
)

LIST ( APPEND ConcurrentCoupler Framework )

CF_ADD_PLUGIN_LIBRARY ( ConcurrentCoupler )

CF_WARN_ORPHAN_FILES()
#ADD_SUBDIRECTORY ( testcases )
