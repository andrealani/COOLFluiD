LIST ( APPEND FiniteVolumePoissonNEQ_files
ComputeBField.cxx
ComputeBField.hh
FiniteVolumePoissonNEQ.hh
PoissonNEQBC.ci
PoissonNEQBC.cxx
PoissonNEQBC.hh
PoissonNEQST.cxx
PoissonNEQST.hh
PoissonCNEQST.cxx
PoissonCNEQST.hh
)

LIST ( APPEND FiniteVolumePoissonNEQ_requires_mods NavierStokes NEQ PoissonNEQ FiniteVolume FiniteVolumeNavierStokes FiniteVolumeNEQ )
LIST ( APPEND FiniteVolumePoissonNEQ_cflibs Framework ShapeFunctions NavierStokes NEQ PoissonNEQ FiniteVolume FiniteVolumeNavierStokes FiniteVolumeNEQ )

CF_ADD_PLUGIN_LIBRARY ( FiniteVolumePoissonNEQ )

CF_WARN_ORPHAN_FILES()
