IF (CF_HAVE_PETSC)

IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "PetscI" )
ELSE ()
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

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

LIST ( APPEND ${MYLIBNAME}_includedirs ${MPI_INCLUDE_DIR} ${PETSC_INCLUDE_DIR} )
LIST ( APPEND ${MYLIBNAME}_libs ${PETSC_LIBRARIES} )
 
IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND ${MYLIBNAME}_cflibs Common Framework )
CF_ADD_PLUGIN_LIBRARY ( PetscI )
ELSE()
 FOREACH (AFILE ${PetscI_files} )
 LIST(APPEND coolfluid-solver_files ../../plugins/Petsc/${AFILE} )
 ENDFOREACH(AFILE)
 CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
 CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_libs  ${coolfluid-solver_libs} )
 
#MESSAGE ("######################### coolfluid-solver_libs in PETSC is ${coolfluid-solver_libs}")

 CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_includedirs  ${coolfluid-solver_includedirs} )
ENDIF()

CF_WARN_ORPHAN_FILES()

ENDIF (CF_HAVE_PETSC)

