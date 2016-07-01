LIST(APPEND FiniteElement_files
FiniteElement.hh
LocalElementData.hh
FElemTypeData.hh
FiniteElementMethod.cxx
FiniteElementMethod.hh
FiniteElementMethodData.cxx
FiniteElementMethodData.hh
ComputeSpaceResidual.cxx
ComputeSpaceResidual.hh
ExplicitComputeSpaceResidual.cxx
ExplicitComputeSpaceResidual.hh
ImplicitComputeSpaceResidual.cxx
ImplicitComputeSpaceResidual.hh
StdComputeTimeResidual.cxx
StdComputeTimeResidual.hh
NewmarkComputeTimeResidual.cxx
NewmarkComputeTimeResidual.hh
ImplicitNewmarkComputeTimeResidual.cxx
ImplicitNewmarkComputeTimeResidual.hh
StdElementComputerStrategy.cxx
StdElementComputerStrategy.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
StdPrepare.cxx
StdPrepare.hh
NumJacobStrategy.cxx
NumJacobStrategy.hh
ComputeJacobStrategy.cxx
ComputeJacobStrategy.hh
FEM_VolumeIntegrator.cxx
FEM_VolumeIntegrator.ci
FEM_VolumeIntegrator.hh
FEM_MeshDataBuilder.cxx
FEM_MeshDataBuilder.hh
NeumannBC.cxx
NeumannBC.hh
DirichletBC.cxx
DirichletBC.hh
NeumannBCImplicitP1Analytical.cxx
NeumannBCImplicitP1Analytical.hh
NeumannBCImplicit.hh
NeumannBCImplicit.cxx
ComputeResidualStrategy.cxx
ComputeResidualStrategy.hh
InitState.cxx
InitState.hh
FEMIntegrableEntity.cxx
FEMIntegrableEntity.hh
ConvectiveEntity.cxx
ConvectiveEntity.hh
DiffusiveEntity.cxx
DiffusiveEntity.hh
NullDiffusiveEntity.cxx
NullDiffusiveEntity.hh
InertiaEntity.cxx
InertiaEntity.hh
NullInertiaEntity.cxx
NullInertiaEntity.hh
IndepSourceEntity.cxx
IndepSourceEntity.hh
LinearSourceEntity.cxx
LinearSourceEntity.hh
NeumannEntity.cxx
NeumannEntity.hh
GalerkinNeumannEntity.cxx
GalerkinNeumannEntity.hh
GalerkinLinearSourceEntity.cxx
GalerkinLinearSourceEntity.hh
GalerkinIndepSourceEntity.cxx
GalerkinIndepSourceEntity.hh
GalerkinConvectiveEntity.cxx
GalerkinConvectiveEntity.hh
ComputeInertiaTerm.cxx
ComputeInertiaTerm.hh
ComputeLinearSourceTerm.cxx
ComputeLinearSourceTerm.hh
ComputeIndepSourceTerm.cxx
ComputeIndepSourceTerm.hh
ComputeDiffusiveTerm.cxx
ComputeDiffusiveTerm.hh
ComputeConvectiveTerm.cxx
ComputeConvectiveTerm.hh
NullComputeConvectiveTerm.cxx
NullComputeConvectiveTerm.hh
NullComputeInertiaTerm.cxx
NullComputeInertiaTerm.hh
NullComputeIndepSourceTerm.cxx
NullComputeIndepSourceTerm.hh
NullComputeLinearSourceTerm.cxx
NullComputeLinearSourceTerm.hh
FEM_HighOrderMeshUpdater.cxx
FEM_HighOrderMeshUpdater.hh
CoupledNeumannBC.hh
CoupledNeumannBC.cxx
CoupledRobinBC.cxx
CoupledRobinBC.hh
CoupledDirichletBC.cxx
CoupledDirichletBC.hh
CoupledNeumannEntity.hh
CoupledNeumannEntity.cxx
GalerkinCoupledNeumannEntity.hh
GalerkinCoupledNeumannEntity.cxx
CoupledRobinImplicitBC.cxx
CoupledRobinImplicitBC.hh
CoupledNeumannImplicitBC.cxx
CoupledNeumannImplicitBC.hh
OutputExtendedConvergence.hh
OutputExtendedConvergence.cxx
MeshMovementDirichletBC.cxx
MeshMovementDirichletBC.hh
)

LIST ( APPEND OPTIONAL_dirfiles
       OutputExtendedConvergenceMPI.hh
       OutputExtendedConvergenceMPI.cxx
     )
IF(CF_HAVE_MPI)
LIST( APPEND FiniteElement_files
      OutputExtendedConvergenceMPI.hh
      OutputExtendedConvergenceMPI.cxx
    )
ENDIF(CF_HAVE_MPI)

LIST ( APPEND FiniteElement_cflibs Framework ShapeFunctions )

CF_ADD_PLUGIN_LIBRARY ( FiniteElement )

########################################################################

LIST(APPEND FiniteElementHeat_files
FiniteElementHeat.hh
ImplicitComputeSpaceResidualP1Analytical.hh
ImplicitComputeSpaceResidualP1Analytical.cxx
GalerkinHeat2DPrimDiffEntity.cxx
GalerkinHeat2DPrimDiffEntity.hh
GalerkinHeat2DPrimCoupledNeumannEntity.cxx
GalerkinHeat2DPrimCoupledNeumannEntity.hh
GalerkinHeat3DPrimDiffEntity.cxx
GalerkinHeat3DPrimDiffEntity.hh
)

LIST ( APPEND FiniteElementHeat_cflibs FiniteElement Heat )
LIST ( APPEND FiniteElementHeat_requires_mods FiniteElement Heat )

CF_ADD_PLUGIN_LIBRARY ( FiniteElementHeat )

########################################################################

LIST(APPEND FiniteElementStructMech_files
FiniteElementStructMech.hh
StructMechComputeStress.hh
StructMechComputeStress.cxx
GalerkinStructMech2DDispDiffEntity.cxx
GalerkinStructMech2DDispDiffEntity.hh
GalerkinStructMech2DAxiDispDiffEntity.cxx
GalerkinStructMech2DAxiDispDiffEntity.hh
GalerkinStructMech3DDispDiffEntity.cxx
GalerkinStructMech3DDispDiffEntity.hh
GalerkinStructMech2DDispInertiaEntity.cxx
GalerkinStructMech2DDispInertiaEntity.hh
GalerkinStructMech3DDispInertiaEntity.cxx
GalerkinStructMech3DDispInertiaEntity.hh
GalerkinStructMech2DDispNeumannEntity.cxx
GalerkinStructMech2DDispNeumannEntity.hh
)

LIST ( APPEND FiniteElementStructMech_cflibs FiniteElement StructMech )
LIST ( APPEND FiniteElementStructMech_requires_mods FiniteElement StructMech )

CF_ADD_PLUGIN_LIBRARY ( FiniteElementStructMech )

########################################################################

LIST(APPEND FiniteElementStructMechHeat_files
FiniteElementStructMechHeat.hh
GalerkinStructMechHeat2DDispDiffEntity.cxx
GalerkinStructMechHeat2DDispDiffEntity.hh
GalerkinStructMechHeat2DAxiDispDiffEntity.cxx
GalerkinStructMechHeat2DAxiDispDiffEntity.hh
GalerkinStructMechHeat2DDispNeumannEntity.cxx
GalerkinStructMechHeat2DDispNeumannEntity.hh
GalerkinStructMechHeat2DDispCoupledNeumannEntity.cxx
GalerkinStructMechHeat2DDispCoupledNeumannEntity.hh
GalerkinStructMechHeat2DIndepSourceEntity.cxx
GalerkinStructMechHeat2DIndepSourceEntity.hh
)

LIST ( APPEND FiniteElementStructMechHeat_cflibs FiniteElementStructMech FiniteElement Heat StructMech StructMechHeat )
LIST ( APPEND FiniteElementStructMechHeat_requires_mods FiniteElementStructMech FiniteElement Heat StructMech StructMechHeat )

CF_ADD_PLUGIN_LIBRARY ( FiniteElementStructMechHeat )

########################################################################

LIST(APPEND FiniteElementChemistry_files
FiniteElementChemistry.hh
GalerkinChemCH43DPrimDiffEntity.cxx
GalerkinChemCH43DPrimDiffEntity.hh
)

LIST ( APPEND FiniteElementChemistry_cflibs FiniteElement CH4 )
LIST ( APPEND FiniteElementChemistry_requires_mods Chemistry )

CF_ADD_PLUGIN_LIBRARY ( FiniteElementChemistry )

########################################################################

CF_WARN_ORPHAN_FILES()

