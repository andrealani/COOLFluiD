IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "Gmsh2CFmesh" )
ELSE () 
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

LIST ( APPEND Gmsh2CFmesh_files
ElementTypeGmsh.hh
Gmsh2CFmesh.hh
Gmsh2CFmeshConverter.cxx
Gmsh2CFmeshConverter.hh
)

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND Gmsh2CFmesh_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( Gmsh2CFmesh )
ELSE()
FOREACH (AFILE ${Gmsh2CFmesh_files} )
LIST(APPEND coolfluid-solver_files ../../plugins/Gmsh2CFmesh/${AFILE} )
ENDFOREACH(AFILE)
CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()
CF_WARN_ORPHAN_FILES()
