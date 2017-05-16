IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "CFmeshExtruder" )
ELSE ()     
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

LIST ( APPEND CFmeshExtruder_files
CFmeshExtruder.hh
Extruder2D.cxx
Extruder2D.hh
Extruder2DFVM.cxx
Extruder2DFVM.hh
Extruder2DDGM.cxx
Extruder2DDGM.hh
)

IF (CF_HAVE_MPI)
LIST ( APPEND CFmeshExtruder_files Extruder2DFVMMPI.cxx Extruder2DFVMMPI.hh )
ENDIF()

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND CFmeshExtruder_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( CFmeshExtruder )
ELSE()
 FOREACH (AFILE ${CFmeshExtruder_files} )
 LIST(APPEND coolfluid-solver_files ../../plugins/CFmeshExtruder/${AFILE} )
 ENDFOREACH(AFILE)
 CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()

CF_WARN_ORPHAN_FILES()
