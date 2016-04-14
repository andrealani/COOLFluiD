LIST ( APPEND NewtonMethod_files
NewtonMethod.hh
SelfAdjustUpdateSol.cxx
SelfAdjustUpdateSol.hh
StdSetup.hh
StdSetup.cxx
StdUnSetup.hh
StdUnSetup.cxx
StdPrepare.hh
StdPrepare.cxx
StdUpdateSol.hh
StdUpdateSol.cxx
TurbUpdateSol.hh
TurbUpdateSol.cxx
NewtonIterator.hh
NewtonIterator.cxx
NewtonIteratorData.hh
NewtonIteratorData.cxx
ResetSystem.hh
ResetSystem.cxx
CopySol.cxx
CopySol.hh
ALE_FVMGeometricAverage.hh
ALE_FVMGeometricAverage.cxx
NewtonIteratorCoupling.hh
NewtonIteratorCoupling.cxx
UpdateSolCoupling.hh
UpdateSolCoupling.cxx
BDF2.hh
BDF2.cxx
BDF2_CN1stStepIntermediate.hh
BDF2_CN1stStepIntermediate.cxx
BDF2_CN1stStepPrepare.hh
BDF2_CN1stStepPrepare.cxx
BDF2_InitCN.hh
BDF2_InitCN.cxx
BDF2Intermediate.hh
BDF2Intermediate.cxx
BDF2Setup.hh
BDF2Setup.cxx
CrankNicholson.hh
CrankNicholson.cxx
CrankNichSetup.hh
CrankNichSetup.cxx
CrankNichUnSetup.cxx
CrankNichUnSetup.hh
CrankNichIntermediate.hh
CrankNichIntermediate.cxx
Linearized.hh
Linearized.cxx
LinearizedBDF2.hh
LinearizedBDF2.cxx
LinearizedBDF2Setup.hh
LinearizedBDF2Setup.cxx
LinearizedSetup.cxx
LinearizedSetup.hh
LinearizedUnSetup.hh
LinearizedUnSetup.cxx
LinearizedPrepare.hh
LinearizedPrepare.cxx
LinearizedIntermediate.hh
LinearizedIntermediate.cxx
LinearizedIntermediateLim.hh
LinearizedIntermediateLim.cxx
NewmarkExplicit.hh
NewmarkExplicit.cxx
NewmarkImplicit.hh
NewmarkImplicit.cxx
NewmarkSetup.hh
NewmarkSetup.cxx
NewmarkPrepare.hh
NewmarkPrepare.cxx
NewmarkResetSystem.cxx
NewmarkResetSystem.hh
NewmarkUnSetup.hh
NewmarkUnSetup.cxx
NewmarkExplicitUpdateSol.hh
NewmarkExplicitUpdateSol.cxx
NewmarkImplicitUpdateSol.hh
NewmarkImplicitUpdateSol.cxx
TwoLayerSetup.hh
TwoLayerSetup.cxx
TwoLayerUnSetup.hh
TwoLayerUnSetup.cxx
TwoLayerPrepare.hh
TwoLayerPrepare.cxx
TwoLayerUpdateSol.hh
TwoLayerUpdateSol.cxx
FSHOSetup.hh
FSHOSetup.cxx
FSHOPrepare.hh
FSHOPrepare.cxx
FSHOUnSetup.hh
FSHOUnSetup.cxx
UFEMSetup.hh
UFEMSetup.cxx
UFEMPrepare.hh
UFEMPrepare.cxx
UFEMUnSetup.hh
UFEMUnSetup.cxx
GReKOUpdateSol.hh
GReKOUpdateSol.cxx
)

LIST ( APPEND OPTIONAL_dirfiles
U.minit.cxx
U.minit.hh
UFEMUpdateSol.cxx
UFEMUpdateSol.hh
)

LIST ( APPEND NewtonMethod_cflibs Framework )


CF_ADD_PLUGIN_LIBRARY ( NewtonMethod )

#######################################################################

LIST ( APPEND NewtonMethodFVM_files
CrankNicholsonLim.hh
CrankNicholsonLim.cxx
CrankNichLimSetup.hh
CrankNichLimSetup.cxx
CrankNichLimUnSetup.cxx
CrankNichLimUnSetup.hh
CrankNichLimInit.hh
CrankNichLimInit.cxx
CrankNichLimIntermediate.hh
CrankNichLimIntermediate.cxx
CrankNichLimPrepare.hh
CrankNichLimPrepare.cxx
UpdateSolFVMCC.cxx
UpdateSolFVMCC.hh
)

LIST ( APPEND NewtonMethodFVM_requires_mods FiniteVolume )

LIST ( APPEND NewtonMethodFVM_cflibs Framework NewtonMethod FiniteVolume )
CF_ADD_PLUGIN_LIBRARY ( NewtonMethodFVM )

#######################################################################

LIST ( APPEND NewtonMethodFVMNavierStokes_files
ImposeHSEquilibriumUpdateSol.cxx
ImposeHSEquilibriumUpdateSol.hh
)

LIST ( APPEND NewtonMethodFVMNavierStokes_requires_mods FiniteVolume NavierStokes)

LIST ( APPEND NewtonMethodFVMNavierStokes_cflibs Framework NewtonMethod FiniteVolume )
CF_ADD_PLUGIN_LIBRARY ( NewtonMethodFVMNavierStokes )

#######################################################################

LIST ( APPEND NewtonMethodMHD_files
NewtonMethodMHD.hh
UpdateSolMHD.hh
UpdateSolMHD.cxx
)

LIST ( APPEND NewtonMethodMHD_requires_mods MHD )
LIST ( APPEND NewtonMethodMHD_cflibs NewtonMethod MHD )


CF_ADD_PLUGIN_LIBRARY ( NewtonMethodMHD )

#######################################################################

CF_WARN_ORPHAN_FILES()
