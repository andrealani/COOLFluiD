IF (CF_HAVE_CGNS)

LIST ( APPEND CGNSWriter_files
CGNSWriter.hh
CGNSWriterAPI.hh
CGWriter.hh
CGWriter.cxx
CGWriterData.hh
CGWriterData.cxx
CGNSHighOrderWriter.cxx
CGNSHighOrderWriter.hh
ParCGNSHighOrderWriter.cxx
ParCGNSHighOrderWriter.hh
StdSetup.hh
StdSetup.cxx
StdUnSetup.hh
StdUnSetup.cxx
)

LIST ( APPEND CGNSWriter_includedirs ${CGNS_INCLUDE_DIR} )
LIST ( APPEND CGNSWriter_libs ${CGNS_LIBRARY}  Framework)

CF_ADD_PLUGIN_LIBRARY ( CGNSWriter )

CF_WARN_ORPHAN_FILES()

ENDIF (CF_HAVE_CGNS)


