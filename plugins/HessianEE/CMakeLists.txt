LIST ( APPEND HessianEE_files
FunctionCalc.cxx
FunctionCalc.hh
HessEE.cxx
HessEE.hh
HessEEData.cxx
HessEEData.hh
HessianEE.hh
IntegralHessCalc.cxx
IntegralHessCalc.hh
StdMetricCalc.cxx
StdMetricCalc.hh
StdSetup.cxx
StdSetup.hh
StdSmoothCom.cxx
StdSmoothCom.hh
StdUnSetup.cxx
StdUnSetup.hh
StdUpdateCom.cxx
StdUpdateCom.hh
)

LIST ( APPEND HessianEE_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( HessianEE )
CF_WARN_ORPHAN_FILES()
