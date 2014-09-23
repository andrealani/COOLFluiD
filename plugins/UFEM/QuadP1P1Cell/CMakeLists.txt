LIST ( APPEND UFEMQuadP1P1Cell_files
UFEMQuadP1P1Cell.hh
CellProps.cxx
CellProps.hh
HeatConduction.cxx
HeatConduction.hh
)

LIST ( APPEND UFEMQuadP1P1Cell_cflibs UFEM )

CF_ADD_PLUGIN_LIBRARY ( UFEMQuadP1P1Cell )

CF_WARN_ORPHAN_FILES()


