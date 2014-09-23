LIST ( APPEND AnalyticalEE_files
AnalyticEE.cxx
AnalyticEE.hh
AnalyticEEData.cxx
AnalyticEEData.hh
AnalyticalEE.hh
ComputeDiscreteError.cxx
ComputeDiscreteError.hh
ComputeIntegralError.cxx
ComputeIntegralError.hh
IntegralErrorFunction.cxx
IntegralErrorFunction.hh
IntegralEntropyErrorEuler.hh
IntegralEntropyErrorEuler.cxx
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
)

LIST ( APPEND AnalyticalEE_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( AnalyticalEE )
CF_WARN_ORPHAN_FILES()
