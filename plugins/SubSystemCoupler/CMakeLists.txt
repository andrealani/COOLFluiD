LIST( APPEND OPTIONAL_dirfiles
  HeatFluxStopCondition.cxx
  HeatFluxStopCondition.hh
)


LIST ( APPEND SubSystemCoupler_files
SubSystemCoupler.hh
Coupler.cxx
Coupler.hh
SubSysCouplerData.cxx
SubSysCouplerData.hh
StdSetup.cxx
StdSetup.hh
FVMCCSetup.cxx
FVMCCSetup.hh
StdPreProcessRead.hh
StdPreProcessRead.cxx
StdPreProcessWrite.cxx
StdPreProcessWrite.hh
FVMCCPreProcessWrite.cxx
FVMCCPreProcessWrite.hh
StructPreProcessWrite.cxx
StructPreProcessWrite.hh
StdMeshMatcherRead.cxx
StdMeshMatcherRead.hh
StdMeshMatcherWrite.cxx
StdMeshMatcherWrite.hh
StdMeshMatcherWrite2.cxx
StdMeshMatcherWrite2.hh
NewtonMeshMatcherWrite.cxx
NewtonMeshMatcherWrite.hh
FVMCCMeshMatcherWrite.cxx
FVMCCMeshMatcherWrite.hh
FVMCCNewtonMeshMatcherWrite.cxx
FVMCCNewtonMeshMatcherWrite.hh
StdReadDataTransfer.cxx
StdReadDataTransfer.hh
StdWriteDataTransfer.cxx
StdWriteDataTransfer.hh
FVMCCReadDataTransfer.cxx
FVMCCReadDataTransfer.hh
FVMCCWriteDataTransfer.cxx
FVMCCWriteDataTransfer.hh
StdPostProcess.cxx
StdPostProcess.hh
PreVariableTransformer.cxx
PreVariableTransformer.hh
PostVariableTransformer.cxx
PostVariableTransformer.hh
NullVariableTransformer.cxx
NullVariableTransformer.hh
NullPreVariableTransformer.cxx
NullPreVariableTransformer.hh
TestVariableTransformer.hh
TestVariableTransformer.cxx
MeshMovementPredictorVariableTransformer.cxx
MeshMovementPredictorVariableTransformer.hh
)

LIST ( APPEND SubSystemCoupler_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( SubSystemCoupler )

##########################################################################

LIST ( APPEND SubSystemCouplerNavierStokes_files
SubSystemCouplerNavierStokes.hh
Euler2DPrimToPressureVariableTransformer.cxx
Euler2DPrimToPressureVariableTransformer.hh
Euler2DConsToPressureVariableTransformer.cxx
Euler2DConsToPressureVariableTransformer.hh
Euler2DConsToPressureVariableTransformerFVMCC.cxx
Euler2DConsToPressureVariableTransformerFVMCC.hh
Euler2DConsToPressureAveragedVariableTransformerFVMCC.cxx
Euler2DConsToPressureAveragedVariableTransformerFVMCC.hh
Euler2DConsToPressurePredictorVariableTransformerFVMCC.cxx
Euler2DConsToPressurePredictorVariableTransformerFVMCC.hh
Euler2DConsToIncompEuler2DdPuvtTransformer.hh
Euler2DConsToIncompEuler2DdPuvtTransformer.cxx
IncompEuler2DdPuvtToTempTransformer.cxx
IncompEuler2DdPuvtToTempTransformer.hh
IncompEuler2DdPuvtToEuler2DConsTransformer.cxx
IncompEuler2DdPuvtToEuler2DConsTransformer.hh
FluidSolidHeatPreVariableTransformerFVMCC.hh
FluidSolidHeatPreVariableTransformerFVMCC.cxx
FluidSolidHeatPreVariableTransformerRhoivtLTEFVMCC.hh
FluidSolidHeatPreVariableTransformerRhoivtLTEFVMCC.cxx
FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC.hh
FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC.cxx
FluidSolidHeatPostVariableTransformerFVMCC.hh
FluidSolidHeatPostVariableTransformerFVMCC.cxx
FluidSolidHeatPreVariableTransformer.hh
FluidSolidHeatPreVariableTransformer.cxx
)

LIST ( APPEND SubSystemCouplerNavierStokes_cflibs SubSystemCoupler NavierStokes )
LIST ( APPEND SubSystemCouplerNavierStokes_requires_mods NavierStokes )

CF_ADD_PLUGIN_LIBRARY ( SubSystemCouplerNavierStokes )

##########################################################################

LIST ( APPEND SubSystemCouplerHeat_files
SubSystemCouplerHeat.hh
Temp2TempAndFluxVariableTransformer.hh
Temp2TempAndFluxVariableTransformer.cxx
TempAndFlux2FluxVariableTransformer.hh
TempAndFlux2FluxVariableTransformer.cxx
TempAndFlux2TempVariableTransformer.hh
TempAndFlux2TempVariableTransformer.cxx
SolidSolidHeatVariableTransformer.hh
SolidSolidHeatVariableTransformer.cxx
ElectroElectroVariableTransformer.hh
ElectroElectroVariableTransformer.cxx
StructMechHeatToDispPreVariableTransformer.hh
StructMechHeatToDispPreVariableTransformer.cxx
StructMechHeatToTempPreVariableTransformer.hh
StructMechHeatToTempPreVariableTransformer.cxx
)

LIST ( APPEND SubSystemCouplerHeat_cflibs SubSystemCoupler Heat )
LIST ( APPEND SubSystemCouplerHeat_requires_mods Heat )

CF_ADD_PLUGIN_LIBRARY ( SubSystemCouplerHeat )


##########################################################################

CF_WARN_ORPHAN_FILES()
#ADD_SUBDIRECTORY ( testcases )
