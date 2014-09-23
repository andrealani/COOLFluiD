IF (CF_HAVE_PETSC)

LIST ( APPEND PetscI_files
BaseSetup.cxx
BaseSetup.hh
BlockJacobiPcJFContext.hh
BlockJacobiPreconditioner.cxx
BlockJacobiPreconditioner.hh
BSORPcContext.hh
BSORPreconditioner.cxx
BSORPreconditioner.hh
DPLURPcJFContext.hh
DPLURPreconditioner.cxx
DPLURPreconditioner.hh
ILUPcContext.hh
ILUPreconditioner.cxx
ILUPreconditioner.hh
JFContext.hh
LUSGSPcJFContext.hh
LUSGSPreconditioner.cxx
LUSGSPreconditioner.hh
MFContext.hh
NewParSetup.cxx
NewParSetup.hh
NullPreconditioner.cxx
NullPreconditioner.hh
ParJFSetup.cxx
ParJFSetup.hh
ParJFSolveSys.cxx
ParJFSolveSys.hh
ParJFSetupGMRESR.cxx
ParJFSetupGMRESR.hh
ParJFSolveSysGMRESR.cxx
ParJFSolveSysGMRESR.hh
ParMFSetup.cxx
ParMFSetup.hh
ParMFSolveSys.cxx
ParMFSolveSys.hh
Petsc.cxx
Petsc.hh
PetscHeaders.hh
PetscIdxMapping.cxx
PetscIdxMapping.hh
PetscLSS.cxx
PetscLSS.hh
PetscLSSData.cxx
PetscLSSData.hh
PetscMatrix.cxx
PetscMatrix.hh
PetscOptions.cxx
PetscOptions.hh
PetscVector.cxx
PetscVector.hh
ShellPreconditioner.cxx
ShellPreconditioner.hh
StdParSolveSys.cxx
StdParSolveSys.hh
StdParUnSetup.cxx
StdParUnSetup.hh
StdSeqSetup.cxx
StdSeqSetup.hh
StdSeqUnSetup.cxx
StdSeqUnSetup.hh
#TriDiagLSS.cxx
#TriDiagLSS.hh
#TriDiagLSSData.cxx
#TriDiagLSSData.hh
TridiagPcJFContext.hh
TridiagPreconditioner.cxx
TridiagPreconditioner.hh
TwoLayerParSetup.cxx
TwoLayerParSetup.hh
TwoLayerParSolveSys.cxx
TwoLayerParSolveSys.hh
TwoLayerSeqSetup.cxx
TwoLayerSeqSetup.hh
TwoLayerSeqSolveSys.cxx
TwoLayerSeqSolveSys.hh
)

LIST ( APPEND PetscI_includedirs ${MPI_INCLUDE_DIR} ${PETSC_INCLUDE_DIR} )

LIST ( APPEND PetscI_cflibs Common Framework ${PETSC_LIBRARIES} )

CF_ADD_PLUGIN_LIBRARY ( PetscI )

CF_WARN_ORPHAN_FILES()

ENDIF (CF_HAVE_PETSC)

