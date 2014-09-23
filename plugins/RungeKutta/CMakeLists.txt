LIST ( APPEND RungeKutta_files
RK.cxx
RK.hh
RKData.cxx
RKData.hh
RungeKutta.hh
RungeKuttaStep.cxx
RungeKuttaStep.hh
StdBackupSol.cxx
StdBackupSol.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
)

LIST ( APPEND RungeKutta_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( RungeKutta )
CF_WARN_ORPHAN_FILES()
