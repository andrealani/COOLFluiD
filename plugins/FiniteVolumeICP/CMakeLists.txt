LIST ( APPEND FiniteVolumeICP_files
FiniteVolumeICP.hh
ICPInductionBC.ci
ICPInductionBC.cxx
ICPInductionBC.hh
ICPInductionEquationSourceTerm.ci
ICPInductionEquationSourceTerm.cxx
ICPInductionEquationSourceTerm.hh
ICPMirror.ci
ICPMirror.cxx
ICPMirror.hh
LorentzForceSourceTermComm.ci
LorentzForceSourceTermComm.cxx
LorentzForceSourceTermComm.hh
LorentzForceSourceTermAlternative.cxx
LorentzForceSourceTermAlternative.hh
RMSJouleHeatSource.cxx
RMSJouleHeatSource.hh
RMSJouleHeatSourceCoupling.ci
RMSJouleHeatSourceCoupling.cxx
RMSJouleHeatSourceCoupling.hh
RMSJouleHeatSourceTerm.cxx
RMSJouleHeatSourceTerm.hh
VectorPotential.cxx
VectorPotential.hh
ICPplasmaFieldComputingBC.cxx
ICPplasmaFieldComputingBC.hh
ICPplasmaFieldComputingBC.ci
StagnationPropsBL.cxx
StagnationPropsBL.hh
)

LIST ( APPEND FiniteVolumeICP_cflibs NavierStokes LTE ICP FiniteVolume FiniteVolumeNavierStokes TecplotWriter )
LIST ( APPEND FiniteVolumeICP_requires_mods NavierStokes LTE ICP FiniteVolume FiniteVolumeNavierStokes TecplotWriter )

LIST ( APPEND FiniteVolumeICPNEQ_files
FiniteVolumeICPNEQ.hh
ICPNEQInductionBC.cxx
ICPNEQInductionEquationSourceTerm.cxx
ICPNEQplasmaFieldComputingBC.cxx
LorentzForceSourceTermCommNEQ.cxx
RMSJouleHeatSourceCouplingNEQ.cxx
)

LIST ( APPEND FiniteVolumeICPNEQ_cflibs  NavierStokes FiniteVolumeICP LTE NEQ ICP FiniteVolume FiniteVolumeNavierStokes FiniteVolumeNEQ TecplotWriter )
LIST ( APPEND FiniteVolumeICPNEQ_requires_mods NavierStokes FiniteVolumeICP LTE  NEQ ICP FiniteVolume FiniteVolumeNavierStokes FiniteVolumeNEQ TecplotWriter )

CF_ADD_PLUGIN_LIBRARY ( FiniteVolumeICP )
CF_ADD_PLUGIN_LIBRARY ( FiniteVolumeICPNEQ )

CF_WARN_ORPHAN_FILES()
