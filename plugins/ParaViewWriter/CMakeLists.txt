IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "ParaViewWriter" )
ELSE ()
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

LIST ( APPEND ParaViewWriter_files
ParaViewWriter.hh
ParaWriter.cxx
ParaWriter.hh
ParaWriterData.cxx
ParaWriterData.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
WriteSolution.cxx
WriteSolution.hh
WriteSolutionHighOrder.cxx
WriteSolutionHighOrder.hh
)

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND ParaViewWriter_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( ParaViewWriter )
ELSE()
 FOREACH (AFILE ${ParaViewWriter_files} )
 LIST(APPEND coolfluid-solver_files ../../plugins/ParaViewWriter/${AFILE} )
 ENDFOREACH(AFILE)
 CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()

CF_WARN_ORPHAN_FILES()
