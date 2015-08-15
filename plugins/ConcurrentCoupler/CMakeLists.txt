LIST ( APPEND ConcurrentCoupler_files
ConcurrentCoupler.hh
ConcurrentCouplerMethod.cxx
ConcurrentCouplerMethod.hh
ConcurrentCouplerData.cxx
ConcurrentCouplerData.hh
StdConcurrentDataTransfer.ci
StdConcurrentDataTransfer.cxx
StdConcurrentDataTransfer.hh
)

LIST ( APPEND ConcurrentCoupler_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( ConcurrentCoupler )

CF_WARN_ORPHAN_FILES()
#ADD_SUBDIRECTORY ( testcases )
