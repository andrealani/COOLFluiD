LIST ( APPEND Gambit2CFmesh_files
ElementTypeGambit.hh
Gambit2CFmesh.hh
Gambit2CFmeshConverter.cxx
Gambit2CFmeshConverter.hh
)

LIST ( APPEND Gambit2CFmesh_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( Gambit2CFmesh )
CF_WARN_ORPHAN_FILES()
