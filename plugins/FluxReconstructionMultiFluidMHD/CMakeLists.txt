LIST( APPEND FluxReconstructionMultiFluidMHD_files
AUSMFluxMultiFluid.cxx
FluxReconstructionMultiFluidMHD.hh
)

LIST ( APPEND FluxReconstructionMultiFluidMHD_cflibs FluxReconstructionMethod MultiFluidMHD Framework)
CF_ADD_PLUGIN_LIBRARY ( FluxReconstructionMultiFluidMHD )

CF_WARN_ORPHAN_FILES()
