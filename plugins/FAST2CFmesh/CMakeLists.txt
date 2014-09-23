LIST ( APPEND FAST2CFmesh_files
FAST2CFmesh.hh
FAST2CFmeshConverter.cxx
FAST2CFmeshConverter.hh
)

LIST ( APPEND FAST2CFmesh_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( FAST2CFmesh )
CF_WARN_ORPHAN_FILES()
