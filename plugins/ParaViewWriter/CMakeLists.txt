LIST ( APPEND ParaViewWriter_files
ParaViewWriter.hh
ParaWriter.cxx
ParaWriter.hh
ParaWriterData.cxx
ParaWriterData.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
WriteSolution.cxx
WriteSolution.hh
WriteSolutionHighOrder.cxx
WriteSolutionHighOrder.hh
)

LIST ( APPEND ParaViewWriter_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( ParaViewWriter )
CF_WARN_ORPHAN_FILES()
