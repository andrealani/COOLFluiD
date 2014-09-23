LIST ( APPEND Dpl2CFmesh_files
Dpl2CFmesh.hh
Dpl2CFmeshConverter.cxx
Dpl2CFmeshConverter.hh
ElementTypeDpl.hh
)

LIST ( APPEND Dpl2CFmesh_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( Dpl2CFmesh )
CF_WARN_ORPHAN_FILES()
