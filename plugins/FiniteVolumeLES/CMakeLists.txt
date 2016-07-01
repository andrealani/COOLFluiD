IF ( CF_HAVE_LESMODELS )

LIST ( APPEND FiniteVolumeLES_files
FiniteVolumeLES.hh
FVMCC_ComputeRHS_LES.hh
FVMCC_ComputeRHS_LES.cxx
FVMCC_ComputeRhsJacob_LES.cxx
FVMCC_ComputeRhsJacob_LES.hh
)

LIST ( APPEND FiniteVolumeLES_cflibs Framework ShapeFunctions FiniteVolume LES)
CF_ADD_PLUGIN_LIBRARY ( FiniteVolumeLES )

CF_WARN_ORPHAN_FILES()

ENDIF ( CF_HAVE_LESMODELS )

