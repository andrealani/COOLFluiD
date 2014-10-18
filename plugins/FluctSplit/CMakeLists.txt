ADD_SUBDIRECTORY (MetaSchemes)

ADD_SUBDIRECTORY (LinEuler)

ADD_SUBDIRECTORY (RANS)

#############################################################################
LIST ( APPEND FluctSplitNavierStokes_files
ArtViscCoeffUnSetup.hh
ArtViscCoeffSetup.cxx
ArtViscCoeffSetup.hh
ArtViscCoeffUnSetup.cxx
WeakSubInletEuler3DConsImpl.hh
WeakSubInletEuler3DConsImpl.cxx
WeakSubInletEuler3DCons.hh
WeakSubInletEuler3DCons.cxx
WeakSubInletEuler3DConsHO.hh
WeakSubInletEuler3DConsHO.cxx
FluctSplitNavierStokes.hh
StrongSubOutletEuler2DCons.hh
StrongSubOutletEuler3DCons.hh
StrongSubOutletEuler3DCons.cxx
WeakSubInletEuler2DCons.hh
WeakFarFieldEuler2DConsAlphaMTP.hh
WeakOutletEuler2DCons.hh
StrongSlipWallEuler3DCons.hh
Euler2DSourceTerm.hh
StrongMirrorEuler2DCons.hh
WeakSlipWallEuler3DCons.hh
WeakSlipWallEuler3DConsHO.hh
StrongMirrorAxisymmEuler2DCons.hh
StrongSlipWallEuler2DCons.hh
StrongSubInletEuler2DCons.hh
StrongSubInletEuler2DConsImpl.hh
StrongSubInletEuler2DConsImpl.cxx
StrongSubInletEuler3DCons.hh
StrongSubInletEuler3DCons.cxx
StrongSubInletEuler3DConsImpl.hh
StrongSubInletEuler3DConsImpl.cxx
WeakSubOutletEuler2DCons.hh
WeakFarFieldEuler2DCons.hh
WeakFarFieldEuler2DConsImpl.hh
WeakSubOutletEuler2DConsImpl.hh
WeakSubInletEuler2DConsImpl.hh
WeakSubInletEuler2DConsPois.hh
WeakSubInletEuler2DConsPois.cxx
StrongMirrorEuler2DConsImpl.hh
WeakSlipWallEuler3DImpl.hh
WeakSubOutletEuler3DConsImpl.hh
WeakSubOutletEuler3DCons.hh
WeakSubOutletEuler3DConsHO.hh
WeakSubOutletEuler3DConsHO.cxx
WeakSlipWallEuler2D.hh
WeakSlipWallEuler2DImpl.hh
StrongNoSlipWallIsothermalNS2DImpl.hh
StrongNoSlipWallIsothermalNS2D.hh
StrongNoSlipWallIsothermalNS3D.hh
StrongNoSlipWallIsothermalNS3DImpl.hh
WeakOutletEuler2DCons.cxx
StrongMirrorEuler2DCons.cxx
StrongSlipWallEuler2DCons.cxx
StrongSubOutletEuler2DCons.cxx
WeakSubInletEuler2DCons.cxx
WeakSlipWallEuler3DConsHO.cxx
WeakSlipWallEuler3DCons.cxx
StrongSubInletEuler2DCons.cxx
StrongSubOutletNonRefSteadyEuler2DCons.cxx
StrongSubOutletNonRefSteadyEuler2DCons.hh
StrongFarFieldNonRefEuler2DCons.cxx
StrongFarFieldNonRefEuler2DCons.hh
StrongSlipWallEuler3DCons.cxx
WeakFarFieldEuler2DConsAlphaMTP.cxx
Euler2DSourceTerm.cxx
StrongMirrorAxisymmEuler2DCons.cxx
WeakSubOutletEuler2DCons.cxx
WeakFarFieldEuler2DCons.cxx
WeakFarFieldEuler2DConsImpl.cxx
WeakSubOutletEuler2DConsImpl.cxx
WeakSubInletEuler2DConsImpl.cxx
StrongMirrorEuler2DConsImpl.cxx
WeakSlipWallEuler3DImpl.cxx
WeakSubOutletEuler3DConsImpl.cxx
WeakSubOutletEuler3DCons.cxx
WeakSlipWallEuler2D.cxx
WeakSlipWallEuler2DImpl.cxx
StrongNoSlipWallIsothermalNS2DImpl.cxx
StrongNoSlipWallIsothermalNS2D.cxx
StrongNoSlipWallIsothermalNS3D.cxx
StrongNoSlipWallIsothermalNS3DImpl.cxx
NavierStokes2DSourceTerm.hh
NavierStokes2DSourceTerm.ci
NavierStokes2DSourceTerm.cxx
NavierStokesTerm.hh
NavierStokesTerm.ci
NavierStokesTerm.cxx
StrongNoSlipWallIsothermalVTImpl.hh
StrongNoSlipWallIsothermalVTImpl.cxx
MixLayerRandomPerturbation.hh
MixLayerRandomPerturbation.cxx
DiffTermCarbuncleFixEuler.hh
DiffTermCarbuncleFixEuler.ci
DiffTermCarbuncleFixEuler.cxx
#DiffTermCarbuncleFixEuler2.hh
#DiffTermCarbuncleFixEuler2.ci
EntropyFix.hh
EntropyFix.cxx
CarbuncleFix.hh
CarbuncleFix.cxx
CarbuncleFix2.hh
CarbuncleFix2.cxx
CarbuncleFixActiveSetup.hh
CarbuncleFixActiveSetup.cxx
CarbuncleFixActiveUnSetup.hh
CarbuncleFixActiveUnSetup.cxx
CbFixDataSetup.cxx
CbFixDataSetup.hh
CbFixDataUnSetup.cxx
CbFixDataUnSetup.hh
BxSchemeSys.cxx
BxSchemeSys.hh
BxSchemeSys.ci
BxSchemeSysLocal.hh
BxSchemeSysLocal.ci
WeakSubInletEuler2DConsRingleb.cxx
WeakSubInletEuler2DConsRingleb.hh
WeakSubOutletEuler2DConsRingleb.cxx
WeakSubOutletEuler2DConsRingleb.hh
UnsteadyWeakSlipWallEuler2DConsImpl.cxx
UnsteadyWeakSlipWallEuler2DConsImpl.hh
UnsteadyWeakSlipWallEuler3DConsImpl.cxx
UnsteadyWeakSlipWallEuler3DConsImpl.hh
UnsteadyWeakSlipWallEuler2DConsHOImpl.cxx
UnsteadyWeakSlipWallEuler2DConsHOImpl.hh
StrongNoSlipWallAdiabaticNS2D.hh
StrongNoSlipWallAdiabaticNS2DImpl.hh
StrongNoSlipWallAdiabaticNS3DImpl.hh
StrongNoSlipWallAdiabaticNS3DImpl.cxx
StrongNoSlipWallAdiabaticNS3D.hh
StrongNoSlipWallAdiabaticNS3D.cxx
StrongNoSlipWallAdiabaticNS2D.cxx
StrongNoSlipWallAdiabaticNS2DImpl.cxx
SaveSourceData.cxx
SaveSourceData.hh
GradComputer.hh
MaximumGradientMagnitude.cxx
MaximumGradientMagnitude.hh
#FluctSplit files
StrongSymmetryAxisImpl.cxx
StrongSymmetryAxisImpl.hh
)

LIST ( APPEND FluctSplitNavierStokes_requires_mods NavierStokes )
LIST ( APPEND FluctSplitNavierStokes_cflibs NavierStokes FluctSplit FluctSplitSystem )

CF_ADD_PLUGIN_LIBRARY ( FluctSplitNavierStokes )

#############################################################################

# LIST ( APPEND FluctSplitAdvectionDiffusion_files
# FluctSplitAdvectionDiffusion.hh
# ScalarDiffusionTerm.hh
# ScalarDiffusionTerm.cxx
# ScalarDiffusionTermHO.hh
# ScalarDiffusionTermHO.cxx
# ScalarLinearadv2DSourceTerm.hh
# ScalarLinearadv2DSourceTerm.cxx
# ScalarLinearadv2DSourceTermHO.hh
# ScalarLinearadv2DSourceTermHO.cxx
# ScalarUnsteady2DSourceTerm.cxx
# ScalarUnsteady2DSourceTerm.hh
# ScalarUnsteady2DSourceTerm2.cxx
# ScalarUnsteady2DSourceTerm2.hh
# ScalarUnsteady2DSourceTermHO.cxx
# ScalarUnsteady2DSourceTermHO.hh
# UnsteadScalarDiffusionTerm.cxx
# UnsteadScalarDiffusionTerm.hh
# LinearAdvDiffHO3D.hh
# LinearAdvDiffHO3D.cxx
# )
# 
# LIST ( APPEND FluctSplitAdvectionDiffusion_requires_mods LinearAdv )
# LIST ( APPEND FluctSplitAdvectionDiffusion_cflibs LinearAdv FluctSplit )
# 
# CF_ADD_PLUGIN_LIBRARY ( FluctSplitAdvectionDiffusion )

#############################################################################

LIST ( APPEND FluctSplitRotationDiffusion_files
FluctSplitRotationDiffusion.hh
ScalarRotationadv2DSourceTerm.hh
ScalarRotationadv2DSourceTerm.cxx
ScalarRotationadv2DSourceTermHO.hh
ScalarRotationadv2DSourceTermHO.cxx
UnsteadScalarRotDiffusionTerm.cxx
UnsteadScalarRotDiffusionTerm.hh
)

LIST ( APPEND FluctSplitRotationDiffusion_requires_mods RotationAdv )
LIST ( APPEND FluctSplitRotationDiffusion_cflibs RotationAdv  FluctSplit )

CF_ADD_PLUGIN_LIBRARY ( FluctSplitRotationDiffusion )

#############################################################################

LIST ( APPEND FluctSplitSpaceTimeNavierStokes_files
FluctSplitSpaceTimeNavierStokes.hh
TwoLayerWeakFarFieldEuler2DConsAlphaMTP.cxx
TwoLayerWeakFarFieldEuler2DConsAlphaMTP.hh
TwoLayerWeakSubInletEuler2DCons.cxx
TwoLayerWeakSubInletEuler2DCons.hh
TwoLayerWeakSubOutletEuler2DCons.cxx
TwoLayerWeakSubOutletEuler2DCons.hh
TwoLayerStrongMirrorEuler2DCons.cxx
TwoLayerStrongMirrorEuler2DCons.hh
TwoLayerWeakSlipWallEuler2DCons.cxx
TwoLayerWeakSlipWallEuler2DCons.hh
StrongSubInletNonRefEuler2DCons.cxx
StrongSubInletNonRefEuler2DCons.hh
StrongSubOutletNonRefEuler2DCons.cxx
StrongSubOutletNonRefEuler2DCons.hh
StrongSubOutletNonRefEuler3DCons.cxx
StrongSubOutletNonRefEuler3DCons.hh
StrongSubOutletNonRefEuler2DConsImpl.cxx
StrongSubOutletNonRefEuler2DConsImpl.hh
StrongSubOutletNonRefEuler3DConsImpl.cxx
StrongSubOutletNonRefEuler3DConsImpl.hh
UnsteadyWeakSubOutletEuler2DConsImpl.cxx
UnsteadyWeakSubOutletEuler2DConsImpl.hh
UnsteadyWeakSubOutletEuler2DCons.cxx
UnsteadyWeakSubOutletEuler2DCons.hh
UnsteadyWeakSlipWallEuler2DImpl.cxx
UnsteadyWeakSlipWallEuler2DImpl.hh
UnsteadyWeakSlipWallEuler2D.cxx
UnsteadyWeakSlipWallEuler2D.hh
UnsteadyWeakSlipWallEulerHO2DCons.cxx
UnsteadyWeakSlipWallEulerHO2DCons.hh
UnsteadyWeakSlipWallEuler3DCons.cxx
UnsteadyWeakSlipWallEuler3DCons.hh
UnsteadyNavierStokesTerm.hh
UnsteadyNavierStokesTerm.cxx
UnsteadyNavierStokesTermquad.hh
UnsteadyNavierStokesTermquad.cxx
UnstP2SUPG_ArtDiffStrategy.hh
UnstP2SUPG_ArtDiffStrategy.cxx
STKT_NlimCSchemeS_nad.hh
STKT_NlimCSchemeS_nad.cxx
RmsSetup.cxx
RmsSetup.hh
QCritSetup.cxx
QCritSetup.hh
DilatationSetup.cxx
DilatationSetup.hh
#VorticitySetup.cxx
#VorticitySetup.hh
RMS_NavierStokes2Dcons.cxx
RMS_NavierStokes2Dcons.hh
RMS_NavierStokes3Dcons.cxx
RMS_NavierStokes3Dcons.hh
QCrit_NavierStokes3Dcons.cxx
QCrit_NavierStokes3Dcons.hh
Dilatation_NavierStokes2Dcons.cxx
Dilatation_NavierStokes2Dcons.hh
Dilatation_NavierStokes3Dcons.cxx
Dilatation_NavierStokes3Dcons.hh
#Vorticity_NavierStokes2Dcons.cxx
#Vorticity_NavierStokes2Dcons.hh
#Vorticity_NavierStokes3Dcons.cxx
#Vorticity_NavierStokes3Dcons.hh
DNS2DSourceTermDampingZone.cxx
DNS2DSourceTermDampingZone.hh
STKT_BDNSSchemeSys.cxx
STKT_BDNSSchemeSys.hh
)

LIST ( APPEND FluctSplitSpaceTimeNavierStokes_requires_mods NavierStokes )
LIST ( APPEND FluctSplitSpaceTimeNavierStokes_cflibs NavierStokes  FluctSplit FluctSplitNavierStokes FluctSplitSpaceTime )

CF_ADD_PLUGIN_LIBRARY ( FluctSplitSpaceTimeNavierStokes )

#############################################################################

LIST ( APPEND FluctSplitMHD_files
FluctSplitMHD.hh
MHD2DProjectionSourceTerm.hh
MHD3DProjectionSourceTerm.hh
StrongSlipWallMHD2D.hh
StrongSlipWallMHD2DProjection.hh
StrongSlipWallMHD3DProjection.hh
StrongSlipWallMHD3D.hh
StrongSlipWallMHD2DImpl.hh
WeakSlipWallMHD2D.hh
WeakSlipWallMHD2DImpl.hh
WeakSlipWallMHD2DProjection.hh
WeakSlipWallMHD2DProjectionImpl.hh
WeakSlipWallMHD3D.hh
WeakSlipWallMHD3DImpl.hh
WeakSlipWallMHD3DProjection.hh
WeakSlipWallMHD3DProjectionImpl.hh
MHD2DSourceTerm.hh
MHD3DSourceTerm.hh
MHD2DSourceTerm.cxx
MHD2DProjectionSourceTerm.cxx
MHD3DProjectionSourceTerm.cxx
WeakSlipWallMHD2D.cxx
WeakSlipWallMHD2DImpl.cxx
WeakSlipWallMHD2DProjection.cxx
WeakSlipWallMHD2DProjectionImpl.cxx
WeakSlipWallMHD3D.cxx
WeakSlipWallMHD3DImpl.cxx
WeakSlipWallMHD3DProjection.cxx
WeakSlipWallMHD3DProjectionImpl.cxx
StrongSlipWallMHD2D.cxx
StrongSlipWallMHD2DProjection.cxx
StrongSlipWallMHD3DProjection.cxx
StrongSlipWallMHD3D.cxx
StrongSlipWallMHD2DImpl.cxx
SuperOutletMHD2DProjection.hh
SuperOutletImplMHD2DProjection.hh
SuperOutletMHD3DProjection.hh
SuperOutletImplMHD3DProjection.hh
SuperInletMHD2DProjection.hh
SuperOutletMHD2DProjection.cxx
SuperOutletImplMHD2DProjection.cxx
SuperOutletMHD3DProjection.cxx
SuperOutletImplMHD3DProjection.cxx
SuperInletMHD2DProjection.cxx
MHD3DSourceTerm.cxx
)

LIST ( APPEND FluctSplitMHD_requires_mods MHD )
LIST ( APPEND FluctSplitMHD_cflibs MHD FluctSplit FluctSplitSystem )

CF_ADD_PLUGIN_LIBRARY ( FluctSplitMHD )

#############################################################################

LIST ( APPEND FluctSplit_files
FluctuationSplit.cxx
FluctuationSplit.hh
FluctuationSplitData.cxx
FluctuationSplitData.hh
InitState.cxx
InitState.hh
StdSetup.cxx
StdSetup.hh
StdSetupIsoP2.cxx
StdSetupIsoP2.hh
StdSetupIsoP3.cxx
StdSetupIsoP3.hh
StdSetupMeshIsoP2CompP1P2.cxx
StdSetupMeshIsoP2CompP1P2.hh
ThetaSetup.cxx
ThetaSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
ThetaUnSetup.cxx
ThetaUnSetup.hh
RDS_MeshDataBuilder.cxx
RDS_MeshDataBuilder.hh
RDS_HighOrderMeshUpdater.cxx
RDS_HighOrderMeshUpdater.hh
NumJacobCouplingStrategy.cxx
NumJacobCouplingStrategy.hh
InwardNormalsData.cxx
InwardNormalsData.hh
ComputeInwardNormals.cxx
ComputeInwardNormals.hh
ComputeInwardNormalsHexaP1.cxx
ComputeInwardNormalsHexaP1.hh
ComputeInwardNormalsPrismP1.cxx
ComputeInwardNormalsPrismP1.hh
ComputeInwardNormalsPyramP1.cxx
ComputeInwardNormalsPyramP1.hh
ComputeInwardNormalsQuadP1.cxx
ComputeInwardNormalsQuadP1.hh
ComputeInwardNormalsTetraP1.cxx
ComputeInwardNormalsTetraP1.hh
ComputeInwardNormalsTriagP1.cxx
ComputeInwardNormalsTriagP1.hh
ComputeInwardNormalsTriagP2.cxx
ComputeInwardNormalsTriagP2.hh
ComputeInwardNormalsTetraP2.cxx
ComputeInwardNormalsTetraP2.hh
ComputeInwardNormalsTriagP3.cxx
ComputeInwardNormalsTriagP3.hh
ComputeRHS.cxx
ComputeRHS.hh
NullComputeSourceTermFSM.cxx
NullComputeSourceTermFSM.hh
NullSplitter.cxx
NullSplitter.hh
CreateBoundaryNodalNormals.cxx
CreateBoundaryNodalNormals.hh
Splitter.cxx
Splitter.hh
FluctuationSplitStrategy.cxx
FluctuationSplitStrategy.hh
RD_SplitStrategy.cxx
RD_SplitStrategy.hh
RD_SplitStrategyMulti.cxx
RD_SplitStrategyMulti.hh
CRD_SplitStrategy.cxx
CRD_SplitStrategy.hh
CRD_SplitStrategyMulti.cxx
CRD_SplitStrategyMulti.hh
ComputeRhsJacob.cxx
ComputeRhsJacob.hh
ComputeRhsJacobCoupling.cxx
ComputeRhsJacobCoupling.hh
ComputeJacobStrategy.cxx
ComputeJacobStrategy.hh
PicardJacobStrategy.cxx
PicardJacobStrategy.hh
NumJacobStrategy.cxx
NumJacobStrategy.hh
NullJacobStrategy.cxx
NullJacobStrategy.hh
StdComputeTimeRhs.cxx
StdComputeTimeRhs.hh
PseudoSteadyTimeRhs.cxx
PseudoSteadyTimeRhs.hh
PseudoSteadyTimeAnalytMatRhs.cxx
PseudoSteadyTimeAnalytMatRhs.hh
NullComputeDiffusiveTerm.cxx
NullComputeDiffusiveTerm.hh
ComputeDiffusiveTerm.cxx
ComputeDiffusiveTerm.hh
P2Normal.cxx
P2Normal.hh
#P2Normal_Tetra.cxx
#P2Normal_Tetra.hh
P3Normal.cxx
P3Normal.hh
ComputeVariablesDerivatives.cxx
ComputeVariablesDerivatives.hh
ArtificialDiffusionStrategy.hh
ArtificialDiffusionStrategy.cxx
Null_ArtDiffStrategy.hh
Null_ArtDiffStrategy.cxx
RDS_SplitterSys.cxx
RDS_SplitterSys.hh
RDS_SplitterSysScalar.cxx
RDS_SplitterSysScalar.hh
RDS_SplitterScalar.cxx
RDS_SplitterScalar.hh
SourceTermSplitter.hh
SourceTermSplitter.cxx
LumpedGalerkinSTSplitter.hh
LumpedGalerkinSTSplitter.cxx
BetaSTSplitter.hh
BetaSTSplitter.cxx
ComputeJacobianFix.hh
ComputeJacobianFix.cxx
NullJacobianFix.hh
NullJacobianFix.cxx
IncludeInFluxSTSplitter.hh
IncludeInFluxSTSplitter.cxx
DistributionData.hh
FluctSplit.hh
ComputeSourceTermFSM.hh
MoveBoundary.hh
MoveBoundary.cxx
SuperInlet.cxx
SuperInlet.hh
SuperInletInterpImpl.cxx
SuperInletInterpImpl.hh
SuperOutlet.cxx
SuperOutlet.hh
WeakBC.cxx
WeakBC.hh
WeakBC2D.cxx
WeakBC2D.hh
WeakBC3D.cxx
WeakBC3D.hh
WeakBC3DHO.cxx
WeakBC3DHO.hh
WeakBC2DImpl.cxx
WeakBC2DImpl.hh
WeakBC3DImpl.cxx
WeakBC3DImpl.hh
WeakSlipWall2DImpl.cxx
WeakSlipWall2DImpl.hh
SuperInletImpl.cxx
SuperInletImpl.hh
CoupledSuperInlet.cxx
CoupledSuperInlet.hh
StrongImplBC.hh
StrongImplBC.cxx
WeakFarField.cxx
WeakFarField.hh
WeakFarField.ci
WeakFarFieldImpl.cxx
WeakFarFieldImpl.hh
WeakFarFieldImpl.ci
WeakSubOutletImpl.cxx
WeakSubOutletImpl.hh
WeakSubOutletImpl.ci
WeakFarFieldInterpImpl.cxx
WeakFarFieldInterpImpl.hh
WeakFarFieldInterpImpl.ci
PeriodicBC.hh
PeriodicBC.cxx
PeriodicBCMPIeqX.cxx
PeriodicBCMPIeqX.hh
PeriodicBC3DMPIeqX.cxx
PeriodicBC3DMPIeqX.hh
PeriodicBC3DMPIeqY.cxx
PeriodicBC3DMPIeqY.hh
PeriodicBCMPIeqY.cxx
PeriodicBCMPIeqY.hh
PeriodicBCImpl.hh
PeriodicBCImpl.cxx
BSysPSICScalar.hh
BSysPSICScalar.cxx
BSchemeCSysScalar.hh
BSchemeCSysScalar.cxx
StrongSymmetryPlane.cxx
StrongSymmetryPlane.hh
StrongSymmetryPlaneImpl.cxx
StrongSymmetryPlaneImpl.hh
)

LIST ( APPEND FluctSplit_cflibs Framework  )

IF ( CF_HAVE_MPI )   
  IF ( CF_HAVE_PARMETIS )
     LIST ( APPEND FluctSplit_files PeriodicBCImplicit.hh PeriodicBCImplicit.cxx ) 
     LIST ( APPEND FluctSplit_includedirs ${PARMETIS_INCLUDE_DIR} )
     LIST ( APPEND FluctSplit_libs ${PARMETIS_LIBRARIES} )
  ENDIF ()
ENDIF()

CF_ADD_PLUGIN_LIBRARY ( FluctSplit )

#############################################################################

LIST ( APPEND FluctSplitSpaceTime_files
FluctSplitSpaceTime.hh
StdALESetup.hh
StdALEUpdate.hh
StdALEPrepare.hh
StdALEUnSetup.hh
StdALESetup.cxx
StdALEUpdate.cxx
StdALEPrepare.cxx
StdALEUnSetup.cxx
TwoLayerComputeRhs.cxx
TwoLayerSetup.cxx
TwoLayerALEUpdate.cxx
TwoLayerUnSetup.cxx
TwoLayerComputeRhs.hh
TwoLayerSetup.hh
TwoLayerALEUpdate.hh
TwoLayerUnSetup.hh
SpaceTimeRD_2LayerSplitStrategy.hh
SpaceTime_Splitter.hh
STKS_NCSchemeScalar.hh
STKS_SplitterSys.hh
STKS_NCSchemeSys.hh
STKS_NlimCSchemeSys.hh
STKS_SplitterScalar.hh
STKS_NlimCSchemeScalar.hh
STKS_LDACSchemeSys.hh
STM_CRD_SplitStrategy.hh
STM_HOCRD_SplitStrategy.hh
STM_HOCRDS2T2_SplitStrategy.hh
STM_RDS_SplitStrategy.hh
STU_CRD_SplitStrategy.hh
STU_QLMJ_SplitStrategy.hh
STU_CRD_SplitStrategyhyb.hh
STU_CRD_SplitStrategyhyb.cxx
STU_CRD_SplitStrategyquad.hh
STU_HOCRD_SplitStrategy.hh
STM_SplitterScalar.hh
STM_SplitterSys.hh
STKT_SplitterScalar.hh
STKT_SplitterSys.hh
STKT_SplitterSyshyb.hh
STKT_SplitterSyshyb.cxx
STM_NSchemeScalar.hh
STM_LDASchemeScalar.hh
STM_LDACSchemeSys.hh
STM_LDASchemeScalarHO.hh
STM_LDASchemeSysHO.hh
STM_LDACSchemeScalar.hh
STKT_LDACSchemeScalar.hh
STKT_LDACSchemeSys.hh
#STKT_LDACSchemeSysScalar.hh
#STKT_LDACSchemeSysScalar.cxx
STKT_NCSchemeSys.hh
STKT_NCSchemeSysquad.hh
STKT_NlimCSchemeSys.hh
STKT_NlimCSchemeScalar.hh
STKT_BCSchemeSys.hh
STKT_BCSchemeSys.cxx
STKT_NCSchemeScalar.hh
STKT_BCSchemeScalar.hh
STM_NlimSchemeScalar.hh
STM_NSchemeSys.hh
STM_LDASchemeSys.hh
STM_NlimSchemeSys.hh
SpaceTime2LayerNSchemeScalar.hh
SpaceTime2LayerNSchemeSys.hh
SpaceTime2LayerNlimSchemeScalar.hh
SpaceTime2LayerNlimSchemeSys.hh
TwoLayerInitState.hh
TwoLayerSuperInlet.hh
TwoLayerWeakBC2D.hh
TwoLayerWeakFarField.hh
TwoLayerComputeTimeRhsImpl.hh
TwoLayerComputeRhsJacobImpl.hh
TwoLayerNumJacobStrategy.hh
TwoLayerSuperInletImpl.hh
SpaceTime_Splitter.cxx
STKS_LDACSchemeSys.cxx
STKS_NCSchemeSys.cxx
STKS_NlimCSchemeScalar.cxx
STKS_SplitterScalar.cxx
STKS_NCSchemeScalar.cxx
STKS_NlimCSchemeSys.cxx
STKS_SplitterSys.cxx
STM_SplitterScalar.cxx
STM_SplitterSys.cxx
STKT_SplitterScalar.cxx
STKT_SplitterScalarQuad.hh
STKT_SplitterScalarQuad.cxx
STKT_SplitterSys.cxx
STKT_SplitterSysQuad.cxx
STKT_SplitterSysQuad.hh
STM_NSchemeScalar.cxx
STM_LDASchemeScalar.cxx
STKT_LDACSchemeScalar.cxx
STKT_LDACSchemeScalarQuad.cxx
STKT_LDACSchemeScalarQuad.hh
STM_LDACSchemeSys.cxx
STM_LDACSchemeScalar.cxx
STKT_LDACSchemeSys.cxx
STKT_LDACSchemeSysQuad.cxx
STKT_LDACSchemeSysQuad.hh
STKT_LDACSchemeSyshyb.cxx
STKT_LDACSchemeSyshyb.hh
STKT_NCSchemeSys.cxx
STKT_NCSchemeSysquad.cxx
STKT_NlimCSchemeSys.cxx
STKT_NCSchemeScalar.cxx
STKT_BCSchemeScalar.cxx
STKT_NlimCSchemeScalar.cxx
STM_NlimSchemeScalar.cxx
STM_NSchemeSys.cxx
STM_BSchemeScalar.hh
STM_BSchemeScalar.cxx
STM_BSchemeSys.hh
STM_BSchemeSys.cxx
STM_LDASchemeSys.cxx
STM_LDASchemeScalarHO.cxx
STM_LDASchemeSysHO.cxx
STM_NlimSchemeSys.cxx
STM_CRD_SplitStrategy.cxx
STM_HOCRD_SplitStrategy.cxx
STM_HOCRDS2T2_SplitStrategy.cxx
STM_RDS_SplitStrategy.cxx
STU_CRD_SplitStrategy.cxx
STU_QLMJ_SplitStrategy.cxx
STU_CRD_SplitStrategyquad.cxx
STU_HOCRD_SplitStrategy.cxx
SpaceTimeRD_2LayerSplitStrategy.cxx
SpaceTime2LayerNSchemeScalar.cxx
SpaceTime2LayerNSchemeSys.cxx
SpaceTime2LayerNlimSchemeScalar.cxx
SpaceTime2LayerNlimSchemeSys.cxx
TwoLayerInitState.cxx
TwoLayerSuperInlet.cxx
TwoLayerWeakBC2D.cxx
TwoLayerWeakFarField.cxx
TwoLayerComputeTimeRhsImpl.cxx
TwoLayerComputeRhsJacobImpl.cxx
TwoLayerNumJacobStrategy.cxx
TwoLayerSuperInletImpl.cxx
STKS_LDACSchemeScalar.cxx
STKS_LDACSchemeScalar.hh
SamplingPoint.hh
SamplingPoint.cxx
)

LIST ( APPEND FluctSplitSpaceTime_cflibs FluctSplit )

CF_ADD_PLUGIN_LIBRARY ( FluctSplitSpaceTime )

#############################################################################

LIST ( APPEND FluctSplitSystem_files
FluctSplitSystem.hh
BSchemeBase.hh
#BLDAFVSchemeCSys.cxx
#BLDAFVSchemeCSys.hh
BLDAHLLSchemeCSys.cxx
BLDAHLLSchemeCSys.hh
#-
#BLWHLLSchemeCSys.cxx
#BLWHLLSchemeCSys.hh
#BLWNSchemeCSys.cxx
#BLWNSchemeCSys.hh
#-
BSchemeCSys.cxx
BSchemeCSys.hh
#BSchemePtSys.cxx
#BSchemePtSys.hh
BSchemeSys.cxx
BSchemeSys.hh
LDASchemeCSys.cxx
LDASchemeCSys.hh
SUPGSchemeCSys_SC.cxx
SUPGSchemeCSys_SC.hh
LDASchemeSys.cxx
LDASchemeSys.hh
NSchemeCSys.cxx
NLimSchemeCSys.cxx
NSchemeCSys.hh
NSchemeCSysScalar.cxx
NSchemeCSysScalar.hh
NLimSchemeCSys.hh
NSchemeSys.cxx
NSchemeSys.hh
NSchemeSysT.ci
NSchemeSysT.cxx
NSchemeSysT.hh
PSISchemeCSys.cxx
PSISchemeCSys.hh
PSISchemeSys.cxx
PSISchemeSys.hh
RDHLLSchemeCSys.cxx
RDHLLSchemeCSys.hh
RusanovSchemeCSys.cxx
RusanovSchemeCSys.hh
BDNSSchemeSys.cxx
BDNSSchemeSys.hh
)

LIST ( APPEND FluctSplitSystem_cflibs FluctSplit )

CF_ADD_PLUGIN_LIBRARY ( FluctSplitSystem )

#############################################################################

LIST ( APPEND FluctSplitScalar_files
FluctSplitScalar.hh
BSchemeCScalar.cxx
BSchemeCScalar.hh
BSchemeScalar.cxx
BSchemeScalar.hh
LDASchemeCScalar.cxx
LDASchemeCScalar.hh
LDASchemeScalar.cxx
LDASchemeScalar.hh
LaxWendroffSchemeCScalar.cxx
LaxWendroffSchemeCScalar.hh
NLimDSchemeCScalar.cxx
NLimDSchemeCScalar.hh
NLimSchemeCScalar.cxx
NLimSchemeCScalar.hh
NSchemeCScalar.cxx
NSchemeCScalar.hh
NSchemeScalar.cxx
NSchemeScalar.hh
#NSchemeScalarT.ci
#NSchemeScalarT.cxx
#NSchemeScalarT.hh
PSISchemeScalar.cxx
PSISchemeScalar.hh
PSISchemeCScalar.cxx
PSISchemeCScalar.hh
RusanovSchemeCScalar.cxx
RusanovSchemeCScalar.hh
SUPGSchemeCScalar.cxx
SUPGSchemeCScalar.hh
#BSchemePtScalar.cxx
#BSchemePtScalar.hh
)

LIST ( APPEND FluctSplitScalar_cflibs FluctSplit )

CF_ADD_PLUGIN_LIBRARY ( FluctSplitScalar )

#############################################################################

#LIST ( APPEND FluctSplitLTE_files
#FluctSplitLTE.hh
#WeakSlipWallEuler2DPrvtLTEImpl.hh
#WeakSlipWallEuler2DPrvtLTEImpl.cxx
#)

#LIST ( APPEND FluctSplitLTE_cflibs FluctSplit LTE )
#LIST ( APPEND FluctSplitLTE_requires_mods FluctSplit LTE FluctSplitSystem )

#CF_ADD_PLUGIN_LIBRARY ( FluctSplitLTE )

#############################################################################

LIST ( APPEND FluctSplitHO_files
FluctSplitHO.hh
HOCRD_SplitStrategy.hh
HOCRD_SplitStrategy.cxx
HOCRD_SplitStrategyIsoP2.hh
HOCRD_SplitStrategyIsoP2.cxx
# HOCRD_SplitStrategyVarBIsoP2.hh
# HOCRD_SplitStrategyVarBIsoP2.cxx
#HOCRD_SplitStrategyIsoP2_GaussLobatto.hh
#HOCRD_SplitStrategyIsoP2_GaussLobatto.cxx
#HOCRD_SplitStrategyIsoP2_GaussLobatto_N4.hh
#HOCRD_SplitStrategyIsoP2_GaussLobatto_N4.cxx
HOCRD_SplitStrategyIsoP3.hh
HOCRD_SplitStrategyIsoP3.cxx
HOCRD_BT_ScalarSplitStrategy.hh
HOCRD_BT_ScalarSplitStrategy.cxx
HOCRDP3_SplitStrategy.hh
HOCRDP3_SplitStrategy.cxx
HOCRD3D.hh
#HOCRD3D_IsoP2.hh
HOCRD3D.cxx
#HOCRD3D_IsoP2.cxx
SUPG_ArtDiffStrategy.hh
SUPG_ArtDiffStrategy.cxx
)



LIST ( APPEND FluctSplitHO_cflibs FluctSplit )
CF_ADD_PLUGIN_LIBRARY ( FluctSplitHO )
#############################################################################

LIST ( APPEND FluctSplitLES_files
FluctSplitLES.hh
UnsteadyGradientTerm.hh
UnsteadyGradientTerm.cxx
UnsteadyClarkTerm.hh
UnsteadyClarkTerm.cxx
UnsteadySmagorinskyTerm.hh
UnsteadySmagorinskyTerm.cxx
UnsteadyWALESTerm.hh
UnsteadyWALESTerm.cxx
UnsteadyLESTerm.hh
UnsteadyLESTerm.cxx
UnsteadyLESTerm.ci
)

LIST ( APPEND FluctSplitLES_cflibs NavierStokes LESvki FluctSplit FluctSplitNavierStokes )
LIST ( APPEND FluctSplitLES_requires_mods NavierStokes LESvki FluctSplit )
CF_ADD_PLUGIN_LIBRARY ( FluctSplitLES )


#############################################################################

ADD_SUBDIRECTORY (HONavierStokes)

#############################################################################

ADD_SUBDIRECTORY (CUDA)

#############################################################################

LIST ( APPEND OPTIONAL_dirfiles 
#  SaveSourceData.cxx
#  LoadSourceData.cxx
#  UpdateCoeffSetup.cxx
#  MixLayer3DRandomPerturbation.cxx
#  UpdateCoeffUnSetup.cxx
#  Statistics_NavierStokes2Dcons.cxx
#  STKT_LDACSchemeSysScalar.cxx
#  STKT_SplitterSysScalar.cxx
#  UpdateCoeffUnSetup.hh
#  SaveSourceData.hh
#  MixLayer3DRandomPerturbation.hh
#  UpdateCoeffSetup.hh
#  LoadSourceData.hh
#  STKT_SplitterSysScalar.hh
#  FluctSplitLinEuler.hh
#  Statistics_NavierStokes2Dcons.hh
#  STKT_LDACSchemeSysScalar.hh
)



CF_WARN_ORPHAN_FILES()
