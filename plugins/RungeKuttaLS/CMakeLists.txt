LIST ( APPEND RungeKuttaLS_files
RKLS.cxx
RKLS.hh
RKLSData.cxx
RKLSData.hh
RungeKuttaLS.hh
RungeKuttaStep.cxx
RungeKuttaStep.hh
StdBackupSol.cxx
StdBackupSol.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
)

LIST ( APPEND RungeKuttaLS_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( RungeKuttaLS )
CF_WARN_ORPHAN_FILES()
