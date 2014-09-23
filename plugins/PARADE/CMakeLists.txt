LIST ( APPEND ParadeI_files
  ParadeLibrary.hh
  ParadeLibrary.cxx
  Parade.hh
)

LIST ( APPEND ParadeI_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( ParadeI )

CF_WARN_ORPHAN_FILES()
