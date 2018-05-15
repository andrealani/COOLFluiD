LIST( APPEND FluxReconstructionNEQ_files
BCFarField.cxx
BCFarField.hh
FluxReconstructionNEQ.hh
CNEQSourceTerm.cxx
CNEQSourceTerm.hh
TNEQSourceTerm.cxx
TNEQSourceTerm.hh
)

LIST ( APPEND FluxReconstructionNEQ_cflibs FluxReconstructionMethod NEQ Framework)
CF_ADD_PLUGIN_LIBRARY ( FluxReconstructionNEQ )

CF_WARN_ORPHAN_FILES()
