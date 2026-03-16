IF (CF_HAVE_PETSC)

LIST ( APPEND FluxReconstructionPetsc_files
FluxReconstructionPetsc.hh
FRP0PcJFContext.hh
FRP0Preconditioner.hh
FRP0Preconditioner.cxx
FRBlockJacobiPcJFContext.hh
FRBlockJacobiPreconditioner.hh
FRBlockJacobiPreconditioner.cxx
FRJacobAssembler.hh
FRJacobAssembler.cxx
)

LIST ( APPEND FluxReconstructionPetsc_requires_mods PetscI FluxReconstructionMethod )
LIST ( APPEND FluxReconstructionPetsc_cflibs Framework PetscI FluxReconstructionMethod )
LIST ( APPEND FluxReconstructionPetsc_includedirs ${MPI_INCLUDE_DIR} ${PETSC_INCLUDE_DIR} )
LIST ( APPEND FluxReconstructionPetsc_libs ${PETSC_LIBRARIES} )

CF_ADD_PLUGIN_LIBRARY ( FluxReconstructionPetsc )

CF_WARN_ORPHAN_FILES()

ENDIF (CF_HAVE_PETSC)
