IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "Tecplot2CFmesh" )
ELSE ()
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

LIST ( APPEND Tecplot2CFmesh_files
ElementTypeTecplot.hh
NodeDim.hh
NodeDim.cxx
Tecplot2CFmesh.hh
Tecplot2CFmeshConverter.cxx
Tecplot2CFmeshConverter.hh
)

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND Tecplot2CFmesh_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( Tecplot2CFmesh )
ELSE()
 FOREACH (AFILE ${Tecplot2CFmesh_files} )
 LIST(APPEND coolfluid-solver_files ../../plugins/Tecplot2CFmesh/${AFILE} )
 ENDFOREACH(AFILE)
 CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()

CF_WARN_ORPHAN_FILES()
