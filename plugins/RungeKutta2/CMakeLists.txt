LIST ( APPEND RungeKutta2_files
CorrectorStep.cxx
CorrectorStep.hh
CorrectorStepTVD.cxx
CorrectorStepTVD.hh
PredictorStep.cxx
PredictorStep.hh
PredictorStepTVD.cxx
PredictorStepTVD.hh
RK2.cxx
RK2.hh
RK2Data.cxx
RK2Data.hh
RungeKutta2.hh
StdBackupSol.cxx
StdBackupSol.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
)

LIST ( APPEND RungeKutta2_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( RungeKutta2 )
CF_WARN_ORPHAN_FILES()
