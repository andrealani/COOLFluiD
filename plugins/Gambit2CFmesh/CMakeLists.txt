IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "Gambit2CFmesh" )
ELSE ()
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

LIST ( APPEND Gambit2CFmesh_files
ElementTypeGambit.hh
Gambit2CFmesh.hh
Gambit2CFmeshConverter.cxx
Gambit2CFmeshConverter.hh
)

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND Gambit2CFmesh_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( Gambit2CFmesh )
ELSE()
FOREACH (AFILE ${Gambit2CFmesh_files} )
LIST(APPEND coolfluid-solver_files ../../plugins/Gambit2CFmesh/${AFILE} )
ENDFOREACH(AFILE)
CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()

CF_WARN_ORPHAN_FILES()
