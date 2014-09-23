LIST ( APPEND TecplotWriterNavierStokes_files
TecplotWriterNavierStokes.hh
WriteInstantAndAvgSolution.cxx
WriteInstantAndAvgSolution.hh
)

LIST ( APPEND TecplotWriterNavierStokes_cflibs TecplotWriter NavierStokes )
LIST ( APPEND TecplotWriterNavierStokes_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( TecplotWriterNavierStokes )
CF_WARN_ORPHAN_FILES()
