LIST ( APPEND LoopMaestro_files
LMaestro.cxx
LMaestro.hh
LoopMaestro.hh
)

LIST ( APPEND LoopMaestro_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( LoopMaestro )
CF_WARN_ORPHAN_FILES()
