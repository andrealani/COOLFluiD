LIST ( APPEND Gmsh2CFmesh_files
ElementTypeGmsh.hh
Gmsh2CFmesh.hh
Gmsh2CFmeshConverter.cxx
Gmsh2CFmeshConverter.hh
)

LIST ( APPEND Gmsh2CFmesh_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( Gmsh2CFmesh )
CF_WARN_ORPHAN_FILES()
