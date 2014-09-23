IF (CF_HAVE_CGNS)

LIST ( APPEND CGNS2CFmesh_files
BaseData.hh
BocoData.hh
CGNS2CFmesh.hh
CGNS2CFmeshConverter.cxx
CGNS2CFmeshConverter.hh
CGNSData.hh
CGNSDefinitions.hh
CGNSException.hh
CGNSReader.cxx
CGNSReader.hh
ElementTypeCGNS.hh
FieldData.hh
GridData.hh
SectionData.hh
SolutionData.hh
ZoneData.hh
)

LIST ( APPEND CGNS2CFmesh_includedirs ${CGNS_INCLUDE_DIR} )
LIST ( APPEND CGNS2CFmesh_libs ${CGNS_LIBRARY}  Framework)

CF_ADD_PLUGIN_LIBRARY ( CGNS2CFmesh )

CF_WARN_ORPHAN_FILES()

ENDIF (CF_HAVE_CGNS)


