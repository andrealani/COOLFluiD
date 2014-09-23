LIST ( APPEND MeshGenerator1D_files
MeshGenerator1D.hh
MeshGenerator1DImpl.cxx
MeshGenerator1DImpl.hh
)

LIST ( APPEND MeshGenerator1D_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( MeshGenerator1D )
CF_WARN_ORPHAN_FILES()
