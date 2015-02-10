LIST ( APPEND MarcoTest_files
DistanceBasedExtrapolator.cxx
GeoDataComputer.cxx
#MarcoTest_BC.cxx
#MarcoTest_BC.hh
MarcoTestAPI.hh
MarcoTest.hh	
MarcoTestMethod.cxx	
MarcoTestMethod.hh	
MarcoTestMethodData.cxx	
MarcoTestMethodData.hh	
Shock.cxx
Shock.hh
StdSetup.cxx
StdUnSetup.cxx
StdUnSetup.hh
)

LIST ( APPEND MarcoTest_cflibs Framework)
CF_ADD_PLUGIN_LIBRARY ( MarcoTest )
CF_WARN_ORPHAN_FILES()
