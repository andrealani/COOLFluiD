LIST ( APPEND THOR2CFmesh_files
CheckNodeNumbering.hh
CheckNodeNumberingHexa.cxx
CheckNodeNumberingHexa.hh
CheckNodeNumberingPrism.cxx
CheckNodeNumberingPrism.hh
CheckNodeNumberingPyram.cxx
CheckNodeNumberingPyram.hh
CheckNodeNumberingTetra.cxx
CheckNodeNumberingTetra.hh
ElementTypeTHOR.hh
THOR2CFmesh.hh
THOR2CFmeshConverter.cxx
THOR2CFmeshConverter.hh
)

LIST ( APPEND THOR2CFmesh_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( THOR2CFmesh )
CF_WARN_ORPHAN_FILES()
