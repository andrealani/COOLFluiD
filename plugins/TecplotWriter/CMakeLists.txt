IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "TecplotWriter" )
ELSE ()
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

LIST ( APPEND TecplotWriter_files
MapGeoEntToTecplot.hh
MapGeoEntToTecplot.cxx
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
TecWriter.cxx
TecWriter.hh
TecWriterData.cxx
TecWriterData.hh
TecplotWriter.hh
TecplotWriterAPI.hh
WriteSolution.cxx
WriteSolution.hh
WriteSolution1D.cxx
WriteSolution1D.hh
WriteSolutionHighOrder.cxx
WriteSolutionHighOrder.hh
WriteSolutionHO.cxx
WriteSolutionHO.hh
WriteSolutionBlock.cxx
WriteSolutionBlock.hh
WriteSolutionBlockDG.cxx
WriteSolutionBlockDG.hh
WriteSolutionBlockFV.cxx
WriteSolutionBlockFV.hh
WriteTecplot.cxx
WriteTecplot.hh
)

IF(CF_HAVE_MPI)
  LIST ( APPEND TecplotWriter_files
    ParWriteSolution.cxx
    ParWriteSolution.hh
    ParWriteSolutionBlock.cxx
    ParWriteSolutionBlock.hh
    )
ENDIF()

LIST ( APPEND OPTIONAL_dirfiles TECXXX.h preutil.c )
IF( CF_HAVE_UNISTD_H )
  LIST ( APPEND TecplotWriter_files TECXXX.h preutil.c )
ENDIF()

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND TecplotWriter_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( TecplotWriter )
ELSE()
FOREACH (AFILE ${TecplotWriter_files} )
LIST(APPEND coolfluid-solver_files ../../plugins/TecplotWriter/${AFILE} )
ENDFOREACH(AFILE)
CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()

CF_WARN_ORPHAN_FILES()
