LIST ( APPEND CFmeshFileWriter_files
CFmeshFileWriter.hh
CFmeshFileWriterAPI.hh
CFmeshWriter.cxx
CFmeshWriter.hh
CFmeshWriterData.cxx
CFmeshWriterData.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
WriteSolution.cxx
WriteSolution.hh
WriteSolutionFluctSplitP2P1.cxx
WriteSolutionFluctSplitP2P1.hh
WriteSolutionDG.cxx
WriteSolutionDG.hh
CFmeshFileWriterDG.ci
CFmeshFileWriterDG.hh
)

IF(CF_HAVE_MPI)
  LIST ( APPEND CFmeshFileWriter_files
ParCFmeshBinaryFileWriter.hh
ParCFmeshBinaryFileWriter.cxx
ParCFmeshFileWriter.hh
ParCFmeshFileWriter.cxx
ParWriteSolution.ci
ParWriteSolution.cxx
ParWriteSolution.hh
)
ENDIF(CF_HAVE_MPI)

LIST ( APPEND CFmeshFileWriter_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( CFmeshFileWriter )

CF_WARN_ORPHAN_FILES()
