LIST ( APPEND Tecplot2CFmesh_files
ElementTypeTecplot.hh
NodeDim.hh
NodeDim.cxx
Tecplot2CFmesh.hh
Tecplot2CFmeshConverter.cxx
Tecplot2CFmeshConverter.hh
)

LIST ( APPEND Tecplot2CFmesh_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( Tecplot2CFmesh )
CF_WARN_ORPHAN_FILES()
