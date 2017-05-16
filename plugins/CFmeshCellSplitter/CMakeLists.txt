IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "CFmeshCellSplitter" )
ELSE ()
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

LIST ( APPEND CFmeshCellSplitter_files
CFmeshCellSplitter.hh
CellSplitter2D.cxx
CellSplitter2D.hh
CellSplitter3D.cxx
CellSplitter3D.hh
CellSplitter2DFVM.cxx
CellSplitter2DFVM.hh
CellSplitter3DFVM.cxx
CellSplitter3DFVM.hh
)

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND CFmeshCellSplitter_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( CFmeshCellSplitter )
ELSE()
 FOREACH (AFILE ${CFmeshCellSplitter_files} )
 LIST(APPEND coolfluid-solver_files ../../plugins/CFmeshCellSplitter/${AFILE} )
 ENDFOREACH(AFILE)
 CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()
