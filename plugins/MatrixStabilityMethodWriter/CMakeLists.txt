LIST ( APPEND MatrixStabilityMethodWriter_files
AddMatrixColumnToFile.cxx
AddMatrixColumnToFile.hh
MatrixStabilityMethod.cxx
MatrixStabilityMethod.hh
MatrixStabilityMethodData.cxx
MatrixStabilityMethodData.hh
MatrixStabilityMethodWriter.hh
SetStates.cxx
SetStates.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
)

LIST ( APPEND MatrixStabilityMethodWriter_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( MatrixStabilityMethodWriter )
CF_WARN_ORPHAN_FILES()
