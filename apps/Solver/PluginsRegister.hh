#include "Environment/FileHandlerInputConcrete.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/DirectFileAccess.hh"
#include "Environment/Environment.hh"
#include "Environment/FileHandlerOutputConcrete.hh"
#include "Environment/DirectFileWrite.hh"

#ifdef CF_HAVE_CURL
#include "Environment/CurlAccessRepository.hh"
#endif

#include "Framework/Framework.hh"
#include "Framework/AbsoluteNormAndMaxIter.hh"
#include "Framework/CellCenteredDiffLocalApproachSparsity.hh"
#include "Framework/CellCenteredSparsity.hh"
#include "Framework/CellVertexSparsity.hh"
#include "Framework/ComputeAllNorms.hh"
#include "Framework/ComputeFaceNormalsHexaP1.hh"
#include "Framework/ComputeFaceNormalsLineP1.hh"
#include "Framework/ComputeFaceNormalsPrismP1.hh"
#include "Framework/ComputeFaceNormalsPyramP1.hh"
#include "Framework/ComputeFaceNormalsQuadP1.hh"
#include "Framework/ComputeFaceNormalsTetraP1.hh"
#include "Framework/ComputeFaceNormalsTriagP1.hh"
#include "Framework/ComputeL2Norm.hh"
#include "Framework/ConstantDTWarn.hh"
#include "Framework/CustomSubSystem.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DetermineCFL.hh"
#include "Framework/ExprComputeCFL.hh"
#include "Framework/ExprComputeDT.hh"
#include "Framework/FVMCC_MeshDataBuilder.hh"
#include "Framework/GlobalMaxNumberStepsCriteria.hh"
#include "Framework/IdentityFilterRHS.hh"
#include "Framework/IdentityFilterState.hh"
#include "Framework/IdentityVarSetMatrixTransformer.hh"
#include "Framework/InteractiveComputeCFL.hh"
#include "Framework/InteractiveComputeDT.hh"
#include "Framework/LMaestro.hh"
#include "Framework/LimitFilterRHS.hh"
#include "Framework/LookupInterpolator.hh"
#include "Framework/MaxComputeDT.hh"
#include "Framework/MaxFilterState.hh"
#include "Framework/MaxNumberStepsCondition.hh"
#include "Framework/MaxNumberSubIterCondition.hh"
#include "Framework/MaxTimeCondition.hh"
#include "Framework/MaxTimeNumberStepsCondition.hh"
#include "Framework/MinMaxFilterState.hh"
#include "Framework/NormAndMaxSubIter.hh"
#include "Framework/NormCondition.hh"
#include "Framework/NullCatalycityModel.hh"
#include "Framework/NullComputeCFL.hh"
#include "Framework/NullComputeDT.hh"
#include "Framework/NullComputeNorm.hh"
#include "Framework/NullContourIntegratorImpl.hh"
#include "Framework/NullConvergenceMethod.hh"
#include "Framework/NullCouplerMethod.hh"
#include "Framework/NullDataProcessing.hh"
#include "Framework/NullDiffusiveVarSet.hh"
#include "Framework/NullDomainModel.hh"
#include "Framework/NullErrorEstimatorMethod.hh"
#include "Framework/NullInertiaVarSet.hh"
#include "Framework/NullIntegrableEntity.hh"
#include "Framework/NullLinearSystemSolver.hh"
#include "Framework/NullLinearizer.hh"
#include "Framework/NullMeshAdapterMethod.hh"
#include "Framework/NullMeshCreator.hh"
#include "Framework/NullMeshDataBuilder.hh"
#include "Framework/NullMeshFormatConverter.hh"
#include "Framework/NullOutputFormatter.hh"
#include "Framework/NullPhysicalModelImpl.hh"
#include "Framework/NullPhysicalPropertyLibrary.hh"
#include "Framework/NullRadiationLibrary.hh"
#include "Framework/NullSourceVarSet.hh"
#include "Framework/NullSpaceMethod.hh"
#include "Framework/NullStateInterpolator.hh"
#include "Framework/NullVarSet.hh"
#include "Framework/NullVarSetMatrixTransformer.hh"
#include "Framework/NullVarSetTransformer.hh"
#include "Framework/NullVolumeIntegratorImpl.hh"
#include "Framework/OnlyMeshSubSystem.hh"
#include "Framework/ParMetis.hh"
#include "Framework/PrePostProcessingSubSystem.hh"
#include "Framework/RelativeNormAndMaxIter.hh"
#include "Framework/RelativeNormAndMaxSubIter.hh"
#include "Framework/SERComputeCFL.hh"
#include "Framework/SMaestro.hh"
#include "Framework/StandardSubSystem.hh"
#include "Framework/SubIterCustomSubSystem.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "CFmeshCellSplitter/CFmeshCellSplitter.hh"
#include "CFmeshCellSplitter/CellSplitter2D.hh"
#include "CFmeshCellSplitter/CellSplitter2DFVM.hh"
#include "CFmeshCellSplitter/CellSplitter3D.hh"
#include "CFmeshCellSplitter/CellSplitter3DFVM.hh"

#include "CFmeshExtruder/CFmeshExtruder.hh"
#include "CFmeshExtruder/Extruder2D.hh"
#include "CFmeshExtruder/Extruder2DDGM.hh"
#include "CFmeshExtruder/Extruder2DFVM.hh"
#include "CFmeshExtruder/Extruder2DFVMMPI.hh"

#include "CFmeshFileReader/CFmeshFileReader.hh"
#include "CFmeshFileReader/CFmeshReader.hh"
#include "CFmeshFileReader/CFmeshReaderData.hh"
#include "CFmeshFileReader/ParReadCFmesh.hh"
#include "CFmeshFileReader/ParCFmeshFileReader.hh"
#include "CFmeshFileReader/ParCFmeshBinaryFileReader.hh"
#include "CFmeshFileReader/ReadCFmesh.hh"
#include "CFmeshFileReader/ReadDummy.hh"
#include "CFmeshFileReader/StdSetup.hh"
#include "CFmeshFileReader/StdUnSetup.hh"

#include "CFmeshFileWriter/CFmeshFileWriter.hh"
#include "CFmeshFileWriter/CFmeshWriter.hh"
#include "CFmeshFileWriter/CFmeshWriterData.hh"
#include "CFmeshFileWriter/ParWriteSolution.hh"
#include "CFmeshFileWriter/ParCFmeshFileWriter.hh"
#include "CFmeshFileWriter/ParCFmeshBinaryFileWriter.hh"
#include "CFmeshFileWriter/StdSetup.hh"
#include "CFmeshFileWriter/StdUnSetup.hh"
#include "CFmeshFileWriter/WriteSolution.hh"
#include "CFmeshFileWriter/WriteSolutionDG.hh"
#include "CFmeshFileWriter/WriteSolutionFluctSplitP2P1.hh"

#include "Gambit2CFmesh/Gambit2CFmesh.hh"
#include "Gambit2CFmesh/Gambit2CFmeshConverter.hh"

#include "Gmsh2CFmesh/Gmsh2CFmeshConverter.hh"
#include "Gmsh2CFmesh/Gmsh2CFmesh.hh"

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/ParWriteSolution.hh"
#include "TecplotWriter/ParWriteSolutionBlock.hh"
#include "TecplotWriter/StdSetup.hh"
#include "TecplotWriter/StdUnSetup.hh"
#include "TecplotWriter/TecWriter.hh"
#include "TecplotWriter/TecWriterData.hh"
#include "TecplotWriter/WriteSolution.hh"
#include "TecplotWriter/WriteSolution1D.hh"
#include "TecplotWriter/WriteSolutionBlock.hh"
#include "TecplotWriter/WriteSolutionBlockDG.hh"
#include "TecplotWriter/WriteSolutionBlockFV.hh"
#include "TecplotWriter/WriteSolutionHO.hh"
#include "TecplotWriter/WriteSolutionHighOrder.hh"

#include "Maxwell/Maxwell.hh"
#include "Maxwell/Maxwell2DCons.hh"
#include "Maxwell/Maxwell2DAdimCons.hh"
#include "Maxwell/Maxwell2DProjectionAdimCons.hh"
#include "Maxwell/Maxwell2DProjectionCons.hh"
#include "Maxwell/Maxwell3DAdimCons.hh"
#include "Maxwell/Maxwell3DCons.hh"
#include "Maxwell/Maxwell3DProjectionCons.hh"
#include "Maxwell/MaxwellModel.hh"
#include "Maxwell/MaxwellModelAdim.hh"
#include "Maxwell/MaxwellProjection.hh"
#include "Maxwell/MaxwellProjectionAdim.hh"

#include "MultiFluidMHD/MultiFluidMHD.hh"
#include "MultiFluidMHD/DiffMFMHD2DHalfRhoiViTi.hh"
#include "MultiFluidMHD/DiffMFMHD2DRhoiViTi.hh"
#include "MultiFluidMHD/DiffMFMHD3DRhoiViTi.hh"
#include "MultiFluidMHD/Euler2DHalfMFMHDConsToRhoiViTiInRhoiViTi.hh"
#include "MultiFluidMHD/Euler2DMFMHDConsToRhoiViTiInRhoiViTi.hh"
#include "MultiFluidMHD/Euler3DMFMHDConsToRhoiViTiInRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD2DCons.hh"
#include "MultiFluidMHD/EulerMFMHD2DConsToRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD2DHalfCons.hh"
#include "MultiFluidMHD/EulerMFMHD2DHalfConsToRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD2DHalfRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD2DHalfRhoiViTiToCons.hh"
#include "MultiFluidMHD/EulerMFMHD2DRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD2DRhoiViTiToCons.hh"
#include "MultiFluidMHD/EulerMFMHD3DCons.hh"
#include "MultiFluidMHD/EulerMFMHD3DConsToRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD3DRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD3DRhoiViTiToCons.hh"
#include "MultiFluidMHD/MultiFluidMHDModel.hh"

#include "ForwardEuler/ForwardEuler.hh"
#include "ForwardEuler/CopySol.hh"
#include "ForwardEuler/FSHOPrepare.hh"
#include "ForwardEuler/FSHOSetup.hh"
#include "ForwardEuler/FSHOUnSetup.hh"
#include "ForwardEuler/FwdEuler.hh"
#include "ForwardEuler/FwdEulerData.hh"
#include "ForwardEuler/StdPrepare.hh"
#include "ForwardEuler/StdSetup.hh"
#include "ForwardEuler/StdUnSetup.hh"
#include "ForwardEuler/TwoLayerPrepare.hh"
#include "ForwardEuler/TwoLayerSetup.hh"
#include "ForwardEuler/TwoLayerUnSetup.hh"
#include "ForwardEuler/TwoLayerUpdateSol.hh"
#include "ForwardEuler/UpdateSol.hh"

#include "NewtonMethod/ALE_FVMGeometricAverage.hh"
#include "NewtonMethod/BDF2.hh"
#include "NewtonMethod/NewtonMethod.hh"
#include "NewtonMethod/NewtonMethodMHD.hh"
#include "NewtonMethod/BDF2Intermediate.hh"
#include "NewtonMethod/BDF2Setup.hh"
#include "NewtonMethod/BDF2_CN1stStepIntermediate.hh"
#include "NewtonMethod/BDF2_CN1stStepPrepare.hh"
#include "NewtonMethod/BDF2_InitCN.hh"
#include "NewtonMethod/CopySol.hh"
#include "NewtonMethod/CrankNichIntermediate.hh"
#include "NewtonMethod/CrankNichLimInit.hh"
#include "NewtonMethod/CrankNichLimIntermediate.hh"
#include "NewtonMethod/CrankNichLimPrepare.hh"
#include "NewtonMethod/CrankNichLimSetup.hh"
#include "NewtonMethod/CrankNichLimUnSetup.hh"
#include "NewtonMethod/CrankNichSetup.hh"
#include "NewtonMethod/CrankNichUnSetup.hh"
#include "NewtonMethod/CrankNicholson.hh"
#include "NewtonMethod/CrankNicholsonLim.hh"
#include "NewtonMethod/FSHOPrepare.hh"
#include "NewtonMethod/FSHOSetup.hh"
#include "NewtonMethod/FSHOUnSetup.hh"
#include "NewtonMethod/GReKOUpdateSol.hh"
#include "NewtonMethod/ImposeHSEquilibriumUpdateSol.hh"
#include "NewtonMethod/Linearized.hh"
#include "NewtonMethod/LinearizedBDF2.hh"
#include "NewtonMethod/LinearizedBDF2Setup.hh"
#include "NewtonMethod/LinearizedIntermediate.hh"
#include "NewtonMethod/LinearizedIntermediateLim.hh"
#include "NewtonMethod/LinearizedPrepare.hh"
#include "NewtonMethod/LinearizedSetup.hh"
#include "NewtonMethod/LinearizedUnSetup.hh"
#include "NewtonMethod/NewmarkExplicit.hh"
#include "NewtonMethod/NewmarkExplicitUpdateSol.hh"
#include "NewtonMethod/NewmarkImplicit.hh"
#include "NewtonMethod/NewmarkImplicitUpdateSol.hh"
#include "NewtonMethod/NewmarkPrepare.hh"
#include "NewtonMethod/NewmarkResetSystem.hh"
#include "NewtonMethod/NewmarkSetup.hh"
#include "NewtonMethod/NewmarkUnSetup.hh"
#include "NewtonMethod/NewtonIterator.hh"
#include "NewtonMethod/NewtonIteratorCoupling.hh"
#include "NewtonMethod/NewtonIteratorData.hh"
#include "NewtonMethod/ResetSystem.hh"
#include "NewtonMethod/SelfAdjustUpdateSol.hh"
#include "NewtonMethod/StdPrepare.hh"
#include "NewtonMethod/StdSetup.hh"
#include "NewtonMethod/StdUnSetup.hh"
#include "NewtonMethod/StdUpdateSol.hh"
#include "NewtonMethod/TurbUpdateSol.hh"
#include "NewtonMethod/TwoLayerPrepare.hh"
#include "NewtonMethod/TwoLayerSetup.hh"
#include "NewtonMethod/TwoLayerUnSetup.hh"
#include "NewtonMethod/TwoLayerUpdateSol.hh"
#include "NewtonMethod/UpdateSolCoupling.hh"
#include "NewtonMethod/UpdateSolFVMCC.hh"
#include "NewtonMethod/UpdateSolMHD.hh"

#include "Petsc/PetscLSSData.hh"
#include "Petsc/Petsc.hh"
#include "Petsc/BSORPreconditioner.hh"
#include "Petsc/BlockJacobiPreconditioner.hh"
#include "Petsc/DPLURPreconditioner.hh"
#include "Petsc/ILUPreconditioner.hh"
#include "Petsc/LUSGSPreconditioner.hh"
#include "Petsc/NewParSetup.hh"
#include "Petsc/NullPreconditioner.hh"
#include "Petsc/ParJFSetup.hh"
#include "Petsc/ParJFSetupGMRESR.hh"
#include "Petsc/ParJFSolveSys.hh"
#include "Petsc/ParJFSolveSysGMRESR.hh"
#include "Petsc/ParMFSolveSys.hh"
#include "Petsc/ParMFSetup.hh"
#include "Petsc/PetscLSS.hh"
#include "Petsc/StdParSolveSys.hh"
#include "Petsc/StdParUnSetup.hh"
#include "Petsc/StdSeqSetup.hh"
#include "Petsc/StdSeqUnSetup.hh"
#include "Petsc/TridiagPreconditioner.hh"
#include "Petsc/TwoLayerParSetup.hh"
#include "Petsc/TwoLayerParSolveSys.hh"
#include "Petsc/TwoLayerSeqSetup.hh"
#include "Petsc/TwoLayerSeqSolveSys.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

class PluginsRegister {
public:

void registerAll() 
{
using namespace Environment;
using namespace Framework;

// Enviroment
#ifdef CF_HAVE_CURL
 Factory<FileHandlerInput>::getInstance().regist
   (new Environment::ObjectProvider<Environment::FileHandlerInputConcrete<  CurlAccessRepository >, 
    Environment::FileHandlerInput,
    Environment::EnvironmentModule >("CurlAccessRepository"));
#endif
 
Factory<FileHandlerInput>::getInstance().regist
(new Environment::ObjectProvider<Environment::FileHandlerInputConcrete< Environment::DirectFileAccess >, 
				 Environment::FileHandlerInput,
				 Environment::EnvironmentModule >("DirectFileAccess"));

Factory<FileHandlerOutput>::getInstance().regist
(new Environment::ObjectProvider<Environment::FileHandlerOutputConcrete< DirectFileWrite >, 
				 Environment::FileHandlerOutput, 
				 Environment::EnvironmentModule >("DirectFileWrite")); 

// Framework
Factory<StopCondition>::getInstance().regist
(new Environment::ObjectProvider<AbsoluteNormAndMaxIter, 
				 StopCondition, FrameworkLib, 1>("AbsoluteNormAndMaxIter"));

Factory<GlobalJacobianSparsity>::getInstance().regist
(new Environment::ObjectProvider<CellCenteredDiffLocalApproachSparsity,
				 GlobalJacobianSparsity, FrameworkLib>("CellCenteredDiffLocalApproach"));

Factory<GlobalJacobianSparsity>::getInstance().regist
(new Environment::ObjectProvider<CellCenteredSparsity,
				 GlobalJacobianSparsity, FrameworkLib>("CellCentered"));

Factory<GlobalJacobianSparsity>::getInstance().regist
(new Environment::ObjectProvider<CellVertexSparsity,
				 GlobalJacobianSparsity, FrameworkLib>("CellVertex"));

Factory<ComputeNorm>::getInstance().regist
(new Environment::ObjectProvider<ComputeAllNorms, ComputeNorm, FrameworkLib, 1>("AllNorms"));

Factory<ComputeNormals>::getInstance().regist
  (new Environment::ObjectProvider<ComputeFaceNormalsHexaP1, ComputeNormals,FrameworkLib>("FaceHexaP1"));
 
 Factory<ComputeNormals>::getInstance().regist
  (new Environment::ObjectProvider<ComputeFaceNormalsLineP1, ComputeNormals,FrameworkLib>("FaceLineP1"));

Factory<ComputeNormals>::getInstance().regist
  (new Environment::ObjectProvider<ComputeFaceNormalsPrismP1, ComputeNormals,FrameworkLib>("FacePrismP1"));

 Factory<ComputeNormals>::getInstance().regist
  (new Environment::ObjectProvider<ComputeFaceNormalsPyramP1, ComputeNormals,FrameworkLib>("FacePyramP1"));
 
Factory<ComputeNormals>::getInstance().regist
  (new Environment::ObjectProvider<ComputeFaceNormalsQuadP1, ComputeNormals,FrameworkLib>("FaceQuadP1"));
 
Factory<ComputeNormals>::getInstance().regist
  (new Environment::ObjectProvider<ComputeFaceNormalsTetraP1, ComputeNormals,FrameworkLib>("FaceTetraP1"));
 
Factory<ComputeNormals>::getInstance().regist
  (new Environment::ObjectProvider<ComputeFaceNormalsTriagP1, ComputeNormals,FrameworkLib>("FaceTriagP1"));

Factory<ComputeNorm>::getInstance().regist
  (new Environment::ObjectProvider<ComputeL2Norm, ComputeNorm, FrameworkLib, 1>("L2"));
 
Factory<ComputeDT>::getInstance().regist
  (new Environment::ObjectProvider<ConstantDTWarn, ComputeDT, FrameworkLib, 1>("ConstantDTWarn"));
 
Factory<SubSystem>::getInstance().regist
  (new Environment::ObjectProvider<CustomSubSystem, SubSystem, FrameworkLib, 1>("CustomSubSystem"));
 
Factory<ComputeNorm>::getInstance().regist
  (new Environment::ObjectProvider<ComputeL2Norm, ComputeNorm, FrameworkLib, 1>("L2"));
 
 Factory<DataProcessingMethod>::getInstance().regist
   (new Environment::ObjectProvider<DataProcessing, DataProcessingMethod, FrameworkLib, 1>("DataProcessing"));
 
 Factory<DataProcessingCom>::getInstance().regist
   (new MethodCommandProvider<NullMethodCommand<DataProcessingData>,
    DataProcessingData, FrameworkLib>("Null"));
 
 Factory<ComputeCFL>::getInstance().regist
   (new Environment::ObjectProvider<DetermineCFL, ComputeCFL, FrameworkLib, 1>("Determine"));

 Factory<ComputeCFL>::getInstance().regist
   (new Environment::ObjectProvider<ExprComputeCFL, ComputeCFL, FrameworkLib, 1>("Function"));

 Factory<ComputeDT>::getInstance().regist
   (new Environment::ObjectProvider<ExprComputeDT, ComputeDT, FrameworkLib, 1>("FunctionDT"));

 Factory<MeshDataBuilder>::getInstance().regist
   (new Environment::ObjectProvider<FVMCC_MeshDataBuilder, MeshDataBuilder, FrameworkLib, 1>("FVMCC"));
 
 Factory<GlobalStopCriteria>::getInstance().regist
   (new Environment::ObjectProvider<GlobalMaxNumberStepsCriteria, GlobalStopCriteria,FrameworkLib, 1>("GlobalMaxNumberSteps"));
 
 Factory<FilterRHS>::getInstance().regist
   (new Environment::ObjectProvider<IdentityFilterRHS,FilterRHS, FrameworkLib,1>("Identity"));
 
 Factory<FilterState>::getInstance().regist
   (new Environment::ObjectProvider<IdentityFilterState,FilterState,FrameworkLib,1>("Identity"));
 
 Factory<VarSetMatrixTransformer>::getInstance().regist
   (new Environment::ObjectProvider<IdentityVarSetMatrixTransformer,VarSetMatrixTransformer,
    FrameworkLib, 1>("Identity"));
 
 Factory<ComputeCFL>::getInstance().regist
   (new Environment::ObjectProvider<InteractiveComputeCFL, ComputeCFL, FrameworkLib, 1>("Interactive"));
 
 Factory<ComputeDT>::getInstance().regist
   (new Environment::ObjectProvider<InteractiveComputeDT, ComputeDT, FrameworkLib, 1>("Interactive"));
 
 Factory<Maestro>::getInstance().regist
   (new ObjectProvider<LMaestro, Maestro, FrameworkLib, 1>("LoopMaestro"));
 
 Factory<FilterRHS>::getInstance().regist
   (new Environment::ObjectProvider<LimitFilterRHS,FilterRHS, FrameworkLib,1>("LimitRHS"));

 Factory<StateInterpolator>::getInstance().regist
   (new ObjectProvider<LookupInterpolator, StateInterpolator, FrameworkLib, 1>("Lookup"));
 
 Factory<ComputeDT>::getInstance().regist
   (new Environment::ObjectProvider<MaxComputeDT, ComputeDT, FrameworkLib, 1>("MaxDT"));

 Factory<FilterState>::getInstance().regist
   (new Environment::ObjectProvider<MaxFilterState,FilterState,FrameworkLib,1>("Max"));
 
 Factory<StopCondition>::getInstance().regist
   (new Environment::ObjectProvider<MaxNumberStepsCondition, 
    StopCondition, FrameworkLib, 1>("MaxNumberSteps"));
 
 Factory<StopCondition>::getInstance().regist
   (new Environment::ObjectProvider<MaxNumberSubIterCondition, 
    StopCondition, FrameworkLib, 1>("MaxNumberSubIter"));
 
 Factory<StopCondition>::getInstance().regist
   (new Environment::ObjectProvider<MaxTimeCondition, 
    StopCondition, FrameworkLib, 1>("MaxTime"));
 
 Factory<StopCondition>::getInstance().regist
   (new Environment::ObjectProvider<MaxTimeNumberStepsCondition, 
    StopCondition, FrameworkLib, 1>("MaxTimeNumberSteps"));
 
 Factory<FilterState>::getInstance().regist
   (new Environment::ObjectProvider<MinMaxFilterState,FilterState,FrameworkLib,1>("MinMax"));
 
 Factory<StopCondition>::getInstance().regist
   (new Environment::ObjectProvider<NormAndMaxSubIter, StopCondition, FrameworkLib, 1>("NormAndMaxSubIter"));
 
 Factory<StopCondition>::getInstance().regist
   (new Environment::ObjectProvider<NormCondition, StopCondition, FrameworkLib, 1>("Norm"));
 
 Factory<CatalycityModel>::getInstance().regist
   (new Environment::ObjectProvider<NullCatalycityModel, CatalycityModel, FrameworkLib, 1>("Null"));
 
 Factory<ComputeCFL>::getInstance().regist
   (new Environment::ObjectProvider<NullComputeCFL, ComputeCFL, FrameworkLib, 1>("Null"));
 
 Factory<ComputeDT>::getInstance().regist
   (new Environment::ObjectProvider<NullComputeDT, ComputeDT, FrameworkLib, 1>("Null"));
 
 Factory<ComputeNorm>::getInstance().regist
   (new Environment::ObjectProvider<NullComputeNorm, ComputeNorm, FrameworkLib, 1>("Null"));

 // IntegratorImplProvider<ContourIntegratorImpl, NullContourIntegratorImpl>()
 
 Factory<ConvergenceMethod>::getInstance().regist
   (new Environment::ObjectProvider<NullConvergenceMethod, ConvergenceMethod, FrameworkLib,1>("Null"));
 
 Factory<CouplerMethod>::getInstance().regist
   (new Environment::ObjectProvider<NullCouplerMethod, CouplerMethod, FrameworkLib,1>("Null"));
 
 Factory<DataProcessingMethod>::getInstance().regist
   (new Environment::ObjectProvider<NullDataProcessing, DataProcessingMethod, FrameworkLib,1>("Null"));

 Factory<DiffusiveVarSet>::getInstance().regist
   (new Environment::ObjectProvider<NullDiffusiveVarSet, DiffusiveVarSet,FrameworkLib,2>("Null"));
 
 Factory<DomainModel>::getInstance().regist
   (new Environment::ObjectProvider<NullDomainModel, DomainModel, FrameworkLib, 1>("Null"));

 Factory<ErrorEstimatorMethod>::getInstance().regist
   (new Environment::ObjectProvider<NullErrorEstimatorMethod, ErrorEstimatorMethod, FrameworkLib,1>("Null"));
 
 Factory<InertiaVarSet>::getInstance().regist
   (new Environment::ObjectProvider<NullInertiaVarSet, InertiaVarSet, FrameworkLib, 1>("Null"));
 
 Factory<IntegrableEntity>::getInstance().regist
   (new Environment::ObjectProvider<NullIntegrableEntity, IntegrableEntity, FrameworkLib>("Null"));
 
 Factory<LinearSystemSolver>::getInstance().regist
   (new Environment::ObjectProvider<NullLinearSystemSolver, LinearSystemSolver, FrameworkLib,1>("Null"));
 
 Factory<JacobianLinearizer>::getInstance().regist
   (new Environment::ObjectProvider<NullLinearizer,JacobianLinearizer,FrameworkLib,1>("Null"));
 
 Factory<MeshAdapterMethod>::getInstance().regist
   (new Environment::ObjectProvider<NullMeshAdapterMethod, MeshAdapterMethod, FrameworkLib,1>("Null"));
 
 Factory<MeshCreator>::getInstance().regist
   (new Environment::ObjectProvider<NullMeshCreator, MeshCreator, FrameworkLib,1>("Null"));
 
 Factory<MeshDataBuilder>::getInstance().regist
   (new Environment::ObjectProvider<NullMeshDataBuilder, MeshDataBuilder, FrameworkLib, 1>("Null"));
 
 Factory<MeshFormatConverter>::getInstance().regist
   (new Environment::ObjectProvider<NullMeshFormatConverter, MeshFormatConverter, FrameworkLib, 1>("Null"));
 
 Factory<OutputFormatter>::getInstance().regist
   (new Environment::ObjectProvider<NullOutputFormatter, OutputFormatter, FrameworkLib,1>("Null"));
 
 Factory<PhysicalModelImpl>::getInstance().regist
   (new Environment::ObjectProvider<NullPhysicalModelImpl, PhysicalModelImpl, FrameworkLib, 1>("Null"));
 
 Factory<PhysicalPropertyLibrary>::getInstance().regist
   (new Environment::ObjectProvider<NullPhysicalPropertyLibrary, PhysicalPropertyLibrary, FrameworkLib, 1>("Null"));
 
 Factory<RadiationLibrary>::getInstance().regist
   (new Environment::ObjectProvider<NullRadiationLibrary, RadiationLibrary, FrameworkLib, 1>("Null"));
 
 Factory<SourceVarSet>::getInstance().regist
   (new Environment::ObjectProvider<NullSourceVarSet,SourceVarSet, FrameworkLib, 1>("Null"));
 
 Factory<SpaceMethod>::getInstance().regist
   (new Environment::ObjectProvider<NullSpaceMethod, SpaceMethod, FrameworkLib,1>("Null"));
 
 Factory<StateInterpolator>::getInstance().regist
   (new ObjectProvider<NullStateInterpolator, StateInterpolator, FrameworkLib, 1>("Null"));
 
 Factory<ConvectiveVarSet>::getInstance().regist
   (new Environment::ObjectProvider<NullVarSet, ConvectiveVarSet, FrameworkLib, 1>("Null"));

 Factory<VarSetMatrixTransformer>::getInstance().regist
   (new Environment::ObjectProvider<NullVarSetMatrixTransformer, VarSetMatrixTransformer, FrameworkLib,1>("Null"));
 
 Factory<VarSetTransformer>::getInstance().regist
   (new Environment::ObjectProvider<NullVarSetTransformer, VarSetTransformer, FrameworkLib,1>("Null"));
 
 // IntegratorImplProvider<VolumeIntegratorImpl, NullVolumeIntegratorImpl>();
 
 Factory<SubSystem>::getInstance().regist
   (new Environment::ObjectProvider<OnlyMeshSubSystem, SubSystem, FrameworkLib, 1>("OnlyMeshSubSystem"));
 
 Factory<MeshPartitioner>::getInstance().regist
   (new Environment::ObjectProvider<ParMetis, MeshPartitioner, FrameworkLib, 1>("ParMetis"));
 
 Factory<SubSystem>::getInstance().regist
   (new Environment::ObjectProvider<PrePostProcessingSubSystem, SubSystem, FrameworkLib, 1>("PrePostProcessingSubSystem"));
 
 Factory<StopCondition>::getInstance().regist
   (new Environment::ObjectProvider<RelativeNormAndMaxIter, StopCondition, FrameworkLib, 1>("RelativeNormAndMaxIter"));
 
 Factory<StopCondition>::getInstance().regist
   (new Environment::ObjectProvider<RelativeNormAndMaxSubIter, StopCondition, FrameworkLib, 1>("RelativeNormAndMaxSubIter"));

 Factory<ComputeCFL>::getInstance().regist
   (new Environment::ObjectProvider<SERComputeCFL, ComputeCFL, FrameworkLib, 1>("SER"));

 Factory<Maestro>::getInstance().regist
   (new ObjectProvider<SMaestro, Maestro, FrameworkLib, 1>("SimpleMaestro"));
 
 Factory<SubSystem>::getInstance().regist
   (new Environment::ObjectProvider<StandardSubSystem, SubSystem, FrameworkLib, 1>("StandardSubSystem"));
 
 Factory<SubSystem>::getInstance().regist
   (new Environment::ObjectProvider<SubIterCustomSubSystem, SubSystem, FrameworkLib, 1>("SubIterCustomSubSystem"));
 
 using namespace IO::CFmeshCellSplitter;
 
 Factory<MeshFormatConverter>::getInstance().regist
   (new Environment::ObjectProvider<CellSplitter2D,  MeshFormatConverter,
    CFmeshCellSplitterModule, 1>("CellSplitter2D"));
 
Factory<MeshFormatConverter>::getInstance().regist
  (new Environment::ObjectProvider<CellSplitter2DFVM,  MeshFormatConverter,
 CFmeshCellSplitterModule, 1>("CellSplitter2DFVM"));
 
 Factory<MeshFormatConverter>::getInstance().regist
   (new Environment::ObjectProvider<CellSplitter3D,  MeshFormatConverter,
    CFmeshCellSplitterModule, 1>("CellSplitter3D"));
 
 Factory<MeshFormatConverter>::getInstance().regist
   (new Environment::ObjectProvider<CellSplitter3DFVM,  MeshFormatConverter,
    CFmeshCellSplitterModule, 1>("CellSplitter3DFVM"));
 
 using namespace IO::CFmeshExtruder;
 
 Factory<MeshFormatConverter>::getInstance().regist
   (new Environment::ObjectProvider<Extruder2D,  MeshFormatConverter,
    CFmeshCellSplitterModule, 1>("Extruder2D"));
 
 Factory<MeshFormatConverter>::getInstance().regist
   (new Environment::ObjectProvider<Extruder2DDGM,  MeshFormatConverter,
    CFmeshCellSplitterModule, 1>("Extruder2DDGM"));
 
 Factory<MeshFormatConverter>::getInstance().regist
   (new Environment::ObjectProvider<Extruder2DFVM,  MeshFormatConverter,
    CFmeshCellSplitterModule, 1>("Extruder2DFVM"));
 
 Factory<MeshFormatConverter>::getInstance().regist
   (new Environment::ObjectProvider<Extruder2DFVMMPI,  MeshFormatConverter,
    CFmeshCellSplitterModule, 1>("Extruder2DFVMMPI"));

 using namespace COOLFluiD::CFmeshFileReader;

Factory<MeshCreator>::getInstance().regist
(new ObjectProvider<CFmeshReader, MeshCreator, CFmeshFileReaderPlugin,1>
 ("CFmeshFileReader"));

Factory<CFmeshReaderCom>::getInstance().regist
(new MethodCommandProvider<NullMethodCommand<CFmeshReaderData>,
 CFmeshReaderData, FrameworkLib>("Null"));

Factory<CFmeshReaderCom>::getInstance().regist
(new MethodCommandProvider<ParReadCFmesh<ParCFmeshFileReader>, 
 CFmeshReaderData, CFmeshFileReaderPlugin>("ParReadCFmesh"));

Factory<CFmeshReaderCom>::getInstance().regist
(new MethodCommandProvider<ParReadCFmesh<ParCFmeshBinaryFileReader>, 
 CFmeshReaderData,CFmeshFileReaderPlugin>("ParReadCFmeshBinary"));

Factory<CFmeshReaderCom>::getInstance().regist
     (new MethodCommandProvider<ReadCFmesh, CFmeshReaderData, 
      CFmeshFileReaderPlugin>("StdReadCFmesh"));

Factory<CFmeshReaderCom>::getInstance().regist
(new MethodCommandProvider<ReadDummy, CFmeshReaderData, CFmeshFileReaderPlugin>("Dummy"));
   
   Factory<CFmeshReaderCom>::getInstance().regist
     (new MethodCommandProvider<COOLFluiD::CFmeshFileReader::StdSetup, CFmeshReaderData, CFmeshFileReaderPlugin>("StdSetup"));
   
   Factory<CFmeshReaderCom>::getInstance().regist
    (new MethodCommandProvider<COOLFluiD::CFmeshFileReader::StdUnSetup, CFmeshReaderData, CFmeshFileReaderPlugin>("StdUnSetup"));  
   
 using namespace CFmeshFileWriter;
 
Factory<OutputFormatter>::getInstance().regist
(new Environment::ObjectProvider<CFmeshWriter, OutputFormatter, CFmeshFileWriterModule,1>
 ("CFmesh"));

Factory<CFmeshWriterCom>::getInstance().regist
(new MethodCommandProvider<NullMethodCommand<CFmeshWriterData>,
 CFmeshWriterData, FrameworkLib>("Null"));

Factory<CFmeshWriterCom>::getInstance().regist
(new MethodCommandProvider<COOLFluiD::CFmeshFileWriter::ParWriteSolution<ParCFmeshFileWriter>, 
 CFmeshWriterData, CFmeshFileWriterModule>("ParWriteSolution"));

Factory<CFmeshWriterCom>::getInstance().regist
  (new MethodCommandProvider<COOLFluiD::CFmeshFileWriter::ParWriteSolution<ParCFmeshBinaryFileWriter>, 
 CFmeshWriterData, CFmeshFileWriterModule>("ParWriteBinarySolution"));
 
Factory<CFmeshWriterCom>::getInstance().regist
  (new MethodCommandProvider<COOLFluiD::CFmeshFileWriter::StdSetup, CFmeshWriterData, CFmeshFileWriterModule>("StdSetup"));

Factory<CFmeshWriterCom>::getInstance().regist
  (new MethodCommandProvider<COOLFluiD::CFmeshFileWriter::StdUnSetup, CFmeshWriterData, CFmeshFileWriterModule>("StdUnSetup"));

 Factory<CFmeshWriterCom>::getInstance().regist
   (new MethodCommandProvider<COOLFluiD::CFmeshFileWriter::WriteSolution, CFmeshWriterData, CFmeshFileWriterModule>("WriteSolution"));

Factory<CFmeshWriterCom>::getInstance().regist
(new MethodCommandProvider<WriteSolutionDG, CFmeshWriterData, CFmeshFileWriterModule>("WriteSolutionDG"));
 
 Factory<CFmeshWriterCom>::getInstance().regist
   (new MethodCommandProvider<WriteSolutionFluctSplitP2P1, CFmeshWriterData, CFmeshFileWriterModule>("WriteSolutionFluctSplitP2P1"));
 
using namespace IO::Gambit2CFmesh;

Factory<MeshFormatConverter>::getInstance().regist
(new Environment::ObjectProvider<Gambit2CFmeshConverter, MeshFormatConverter, 
 Gambit2CFmeshModule, 1>("Gambit2CFmesh"));

 using namespace IO::Gmsh2CFmesh;
 
 Factory<MeshFormatConverter>::getInstance().regist
(new Environment::ObjectProvider<Gmsh2CFmeshConverter, MeshFormatConverter, 
 Gmsh2CFmeshModule, 1>("Gmsh2CFmesh"));
 
 using namespace TecplotWriter;
 
 Factory<TecWriterCom>::getInstance().regist
   (new MethodCommandProvider<COOLFluiD::TecplotWriter::ParWriteSolution, TecWriterData, TecplotWriterModule>
    ("ParWriteSolution"));
 
 Factory<TecWriterCom>::getInstance().regist
   (new MethodCommandProvider<ParWriteSolutionBlock, TecWriterData, TecplotWriterModule>
    ("ParWriteSolutionBlock"));
 
 Factory<TecWriterCom>::getInstance().regist
   (new MethodCommandProvider<COOLFluiD::TecplotWriter::StdSetup, TecWriterData, TecplotWriterModule>
    ("StdSetup"));
 
 Factory<TecWriterCom>::getInstance().regist
   (new MethodCommandProvider<COOLFluiD::TecplotWriter::StdUnSetup, TecWriterData, TecplotWriterModule>
    ("StdUnSetup"));
 
 Factory<OutputFormatter>::getInstance().regist
   (new Environment::ObjectProvider<TecWriter, OutputFormatter, TecplotWriterModule,1>
    ("Tecplot"));
 
 Factory<TecWriterCom>::getInstance().regist
   (new MethodCommandProvider<NullMethodCommand<TecWriterData>, TecWriterData, 
    TecplotWriterModule>("Null"));
 
 Factory<TecWriterCom>::getInstance().regist
   (new MethodCommandProvider<COOLFluiD::TecplotWriter::WriteSolution, TecWriterData, TecplotWriterModule>
    ("WriteSolution"));
 
 Factory<TecWriterCom>::getInstance().regist
   (new MethodCommandProvider<WriteSolution1D, TecWriterData, TecplotWriterModule>
    ("WriteSolution1D"));

 Factory<TecWriterCom>::getInstance().regist
   (new MethodCommandProvider<WriteSolutionBlock, TecWriterData, TecplotWriterModule>
    ("WriteSolutionBlock"));

 Factory<TecWriterCom>::getInstance().regist
   (new MethodCommandProvider<WriteSolutionBlockDG, TecWriterData, TecplotWriterModule>
    ("WriteSolutionBlockDG"));
 
 Factory<TecWriterCom>::getInstance().regist
   (new MethodCommandProvider<WriteSolutionBlockFV, TecWriterData, TecplotWriterModule>
    ("WriteSolutionBlockFV"));
 
 Factory<TecWriterCom>::getInstance().regist
   (new MethodCommandProvider<WriteSolutionHO, TecWriterData, TecplotWriterModule>
    ("WriteSolutionHO"));
 
 Factory<TecWriterCom>::getInstance().regist
   (new MethodCommandProvider<WriteSolutionHighOrder, TecWriterData, TecplotWriterModule>
    ("WriteSolutionHighOrder")); 
 
 using namespace Physics::Maxwell;
 
 Factory<ConvectiveVarSet>::getInstance().regist
   (new Environment::ObjectProvider<Maxwell2DAdimCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell2DAdimCons"));
 
Factory<ConvectiveVarSet>::getInstance().regist
   (new Environment::ObjectProvider<Maxwell2DCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell2DCons"));
 
Factory<ConvectiveVarSet>::getInstance().regist
(new Environment::ObjectProvider<Maxwell2DProjectionAdimCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell2DProjectionAdimCons"));

Factory<ConvectiveVarSet>::getInstance().regist
(new Environment::ObjectProvider<Maxwell2DProjectionCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell2DProjectionCons"));

Factory<ConvectiveVarSet>::getInstance().regist
(new Environment::ObjectProvider<Maxwell3DAdimCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell3DAdimCons"));

Factory<ConvectiveVarSet>::getInstance().regist
(new Environment::ObjectProvider<Maxwell3DCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell3DCons"));

Factory<ConvectiveVarSet>::getInstance().regist
(new Environment::ObjectProvider<Maxwell3DProjectionCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell3DProjectionCons"));

Factory<PhysicalModelImpl>::getInstance().regist
(new Environment::ObjectProvider<MaxwellModel<DIM_2D>, PhysicalModelImpl,MaxwellModule, 1>("Maxwell2D"));

Factory<PhysicalModelImpl>::getInstance().regist
(new Environment::ObjectProvider<MaxwellModel<DIM_3D>, PhysicalModelImpl,MaxwellModule, 1>("Maxwell3D"));

Factory<PhysicalModelImpl>::getInstance().regist
(new Environment::ObjectProvider<MaxwellModelAdim<DIM_2D>, PhysicalModelImpl, MaxwellModule,1>("Maxwell2DAdim"));

Factory<PhysicalModelImpl>::getInstance().regist
(new Environment::ObjectProvider<MaxwellModelAdim<DIM_3D>, PhysicalModelImpl, MaxwellModule,1>("Maxwell3DAdim"));

Factory<PhysicalModelImpl>::getInstance().regist
(new Environment::ObjectProvider<MaxwellProjection<DIM_2D>, PhysicalModelImpl, MaxwellModule,1>("Maxwell2DProjection"));

Factory<PhysicalModelImpl>::getInstance().regist
(new Environment::ObjectProvider<MaxwellProjection<DIM_3D>, PhysicalModelImpl, MaxwellModule,1>("Maxwell3DProjection"));

Factory<PhysicalModelImpl>::getInstance().regist
(new Environment::ObjectProvider<MaxwellProjectionAdim<DIM_2D>, PhysicalModelImpl, MaxwellModule,1>("Maxwell2DProjectionAdim"));

Factory<PhysicalModelImpl>::getInstance().regist
(new Environment::ObjectProvider<MaxwellProjectionAdim<DIM_3D>, PhysicalModelImpl, MaxwellModule,1>("Maxwell3DProjectionAdim"));
 
 using namespace Physics::MultiFluidMHD;
 
 Factory<DiffusiveVarSet>::getInstance().regist
  (new Environment::ObjectProvider<DiffMFMHD2DHalfRhoiViTi, DiffusiveVarSet, MultiFluidMHDModule, 2>("MultiFluidMHD2DHalfRhoiViTi"));
 
Factory<DiffusiveVarSet>::getInstance().regist
(new Environment::ObjectProvider<DiffMFMHD2DRhoiViTi, DiffusiveVarSet, MultiFluidMHDModule, 2>("MultiFluidMHD2DRhoiViTi"));
 
Factory<DiffusiveVarSet>::getInstance().regist
(new Environment::ObjectProvider<DiffMFMHD3DRhoiViTi, DiffusiveVarSet, MultiFluidMHDModule, 2>("MultiFluidMHD3DRhoiViTi"));

 Factory<VarSetMatrixTransformer>::getInstance().regist
   (new Environment::ObjectProvider<Euler2DHalfMFMHDConsToRhoiViTiInRhoiViTi, VarSetMatrixTransformer, 
    MultiFluidMHDModule, 1> ("Euler2DHalfMFMHDConsToRhoiViTiInRhoiViTi"));
 
 Factory<VarSetMatrixTransformer>::getInstance().regist
   (new Environment::ObjectProvider<Euler2DMFMHDConsToRhoiViTiInRhoiViTi, VarSetMatrixTransformer, 
    MultiFluidMHDModule, 1> ("Euler2DMFMHDConsToRhoiViTiInRhoiViTi"));
 
 Factory<VarSetMatrixTransformer>::getInstance().regist
   (new Environment::ObjectProvider<Euler3DMFMHDConsToRhoiViTiInRhoiViTi, VarSetMatrixTransformer, 
    MultiFluidMHDModule, 1>("Euler3DMFMHDConsToRhoiViTiInRhoiViTi"));
 
 Factory<ConvectiveVarSet>::getInstance().regist
   (new Environment::ObjectProvider<EulerMFMHD2DCons, ConvectiveVarSet, MultiFluidMHDModule, 1>
    ("EulerMFMHD2DCons"));
 
 Factory<VarSetTransformer>::getInstance().regist
   (new Environment::ObjectProvider<EulerMFMHD2DConsToRhoiViTi, VarSetTransformer, MultiFluidMHDModule, 1>
    ("EulerMFMHD2DConsToRhoiViTi"));
 
 Factory<ConvectiveVarSet>::getInstance().regist
   (new Environment::ObjectProvider<EulerMFMHD2DHalfCons, ConvectiveVarSet, MultiFluidMHDModule, 1>
    ("EulerMFMHD2DHalfCons"));
 
 Factory<VarSetTransformer>::getInstance().regist
   (new Environment::ObjectProvider<EulerMFMHD2DHalfConsToRhoiViTi, VarSetTransformer, MultiFluidMHDModule, 1>
    ("EulerMFMHD2DHalfConsToRhoiViTi"));
 
 Factory<ConvectiveVarSet>::getInstance().regist
   (new Environment::ObjectProvider<EulerMFMHD2DHalfRhoiViTi, ConvectiveVarSet, MultiFluidMHDModule, 1>
    ("EulerMFMHD2DHalfRhoiViTi"));
 
Factory<VarSetTransformer>::getInstance().regist
(new Environment::ObjectProvider<EulerMFMHD2DHalfRhoiViTiToCons, VarSetTransformer, MultiFluidMHDModule, 1>
 ("EulerMFMHD2DHalfRhoiViTiToCons"));

Factory<ConvectiveVarSet>::getInstance().regist
(new Environment::ObjectProvider<EulerMFMHD2DRhoiViTi, ConvectiveVarSet, MultiFluidMHDModule, 1>
 ("EulerMFMHD2DRhoiViTi"));

Factory<VarSetTransformer>::getInstance().regist
(new Environment::ObjectProvider<EulerMFMHD2DRhoiViTiToCons, VarSetTransformer, MultiFluidMHDModule, 1>
 ("EulerMFMHD2DRhoiViTiToCons"));

Factory<ConvectiveVarSet>::getInstance().regist
(new Environment::ObjectProvider<EulerMFMHD3DCons, ConvectiveVarSet, MultiFluidMHDModule, 1>
 ("EulerMFMHD3DCons"));

Factory<VarSetTransformer>::getInstance().regist
  (new Environment::ObjectProvider<EulerMFMHD3DConsToRhoiViTi, VarSetTransformer, MultiFluidMHDModule, 1>
 ("EulerMFMHD3DConsToRhoiViTi"));
 
 Factory<ConvectiveVarSet>::getInstance().regist
   (new Environment::ObjectProvider<EulerMFMHD3DRhoiViTi, ConvectiveVarSet, MultiFluidMHDModule, 1>
    ("EulerMFMHD3DRhoiViTi"));
 
 Factory<VarSetTransformer>::getInstance().regist
   (new Environment::ObjectProvider<EulerMFMHD3DRhoiViTiToCons, VarSetTransformer, MultiFluidMHDModule, 1>
    ("EulerMFMHD3DRhoiViTiToCons"));
 
 Factory<PhysicalModelImpl>::getInstance().regist
   (new Environment::ObjectProvider<MultiFluidMHDModel<DIM_2D>, PhysicalModelImpl,MultiFluidMHDModule, 1>
    ("MultiFluidMHD2D"));
 
 Factory<PhysicalModelImpl>::getInstance().regist
   (new Environment::ObjectProvider<MultiFluidMHDModel<DIM_3D>, PhysicalModelImpl,MultiFluidMHDModule, 1>
    ("MultiFluidMHD3D"));
 
 using namespace ForwardEuler;
 
 Factory<FwdEulerCom>::getInstance().regist
   (new MethodCommandProvider<COOLFluiD::ForwardEuler::CopySol,FwdEulerData,ForwardEulerLib >("StdCopySol"));
 
 Factory<FwdEulerCom>::getInstance().regist
   (new MethodCommandProvider<FSHOPrepare, FwdEulerData, ForwardEulerLib>("FSHOPrepare"));
 
 Factory<FwdEulerCom>::getInstance().regist
   (new MethodCommandProvider<FSHOSetup, FwdEulerData, ForwardEulerLib>("FSHOSetup"));
 
 Factory<FwdEulerCom>::getInstance().regist
   (new MethodCommandProvider<FSHOUnSetup, FwdEulerData, ForwardEulerLib>("FSHOUnSetup"));
 
 Factory<ConvergenceMethod>::getInstance().regist
   (new Environment::ObjectProvider<FwdEuler, ConvergenceMethod, ForwardEulerLib, 1>
    ("FwdEuler"));

Factory<FwdEulerCom>::getInstance().regist
(new MethodCommandProvider<NullMethodCommand<FwdEulerData>, FwdEulerData, ForwardEulerLib>
 ("Null"));

 Factory<FwdEulerCom>::getInstance().regist
   (new MethodCommandProvider<StdPrepare, FwdEulerData, ForwardEulerLib>
    ("StdPrepare"));
 
 Factory<FwdEulerCom>::getInstance().regist
   (new MethodCommandProvider<COOLFluiD::ForwardEuler::StdSetup, FwdEulerData, ForwardEulerLib>("StdSetup"));

Factory<FwdEulerCom>::getInstance().regist
(new MethodCommandProvider<COOLFluiD::ForwardEuler::StdUnSetup, FwdEulerData, ForwardEulerLib>("StdUnSetup"));

Factory<FwdEulerCom>::getInstance().regist
(new MethodCommandProvider<TwoLayerPrepare, FwdEulerData, ForwardEulerLib>
 ("TwoLayerPrepare"));

Factory<FwdEulerCom>::getInstance().regist
(new MethodCommandProvider<TwoLayerSetup, FwdEulerData, ForwardEulerLib>("TwoLayerSetup"));

Factory<FwdEulerCom>::getInstance().regist
(new MethodCommandProvider<TwoLayerUnSetup, FwdEulerData, ForwardEulerLib>("TwoLayerUnSetup"));

Factory<FwdEulerCom>::getInstance().regist
(new MethodCommandProvider<TwoLayerUpdateSol, FwdEulerData, ForwardEulerLib>
 ("TwoLayerUpdateSol"));

Factory<FwdEulerCom>::getInstance().regist
  (new MethodCommandProvider<COOLFluiD::ForwardEuler::UpdateSol, FwdEulerData, ForwardEulerLib>
 ("StdUpdateSol"));
 
using namespace Numerics::NewtonMethod;
 
Factory<NewtonIteratorCom>::getInstance().regist
  (new MethodCommandProvider<ALE_FVMGeometricAverage, NewtonIteratorData, NewtonMethodModule> ("ALE_FVMGeometricAverage"));

 Factory<ConvergenceMethod>::getInstance().regist
   (new Environment::ObjectProvider<BDF2, ConvergenceMethod, NewtonMethodModule, 1>
    ("BDF2"));

 Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<BDF2Intermediate, NewtonIteratorData, NewtonMethodModule>
    ("BDF2Intermediate"));

Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<BDF2Setup, NewtonIteratorData, NewtonMethodModule>
    ("BDF2Setup"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<BDF2_CN1stStepIntermediate, NewtonIteratorData, NewtonMethodModule>
 ("CN1stStepIntermediate"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<BDF2_CN1stStepPrepare, NewtonIteratorData, NewtonMethodModule>
 ("CN1stStepPrepare"));

Factory<ConvergenceMethod>::getInstance().regist
   (new Environment::ObjectProvider<BDF2_InitCN, ConvergenceMethod, NewtonMethodModule, 1>
    ("BDF2_InitCN"));

 Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::CopySol, NewtonIteratorData, NewtonMethodModule>("CopySol"));
 
 Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<CrankNichIntermediate, NewtonIteratorData, NewtonMethodModule> 
    ("CrankNichIntermediate"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<CrankNichLimInit, NewtonIteratorData, NewtonMethodModule>
 ("CrankNichLimInit"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<CrankNichLimIntermediate, NewtonIteratorData, NewtonMethodModule>
 ("CrankNichLimIntermediate"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<CrankNichLimPrepare, NewtonIteratorData, NewtonMethodModule>
 ("CrankNichLimPrepare"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<CrankNichLimSetup, NewtonIteratorData, NewtonMethodModule>
 ("CrankNichLimSetup"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<CrankNichLimUnSetup, NewtonIteratorData, NewtonMethodModule> 
 ("CrankNichLimUnSetup"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<CrankNichSetup, NewtonIteratorData, NewtonMethodModule>
 ("CrankNichSetup"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<CrankNichUnSetup, NewtonIteratorData, NewtonMethodModule>
 ("CrankNichUnSetup"));
 
Factory<ConvergenceMethod>::getInstance().regist
   (new Environment::ObjectProvider<CrankNicholson, ConvergenceMethod, NewtonMethodModule, 1>
    ("CrankNicholson"));

Factory<ConvergenceMethod>::getInstance().regist
   (new Environment::ObjectProvider<CrankNicholsonLim, ConvergenceMethod, NewtonMethodModule, 1>
   ("CrankNicholsonLim"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::FSHOPrepare, NewtonIteratorData, NewtonMethodModule>
("FSHOPrepare"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::FSHOSetup, NewtonIteratorData, NewtonMethodModule>
("FSHOSetup"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::FSHOUnSetup, NewtonIteratorData, NewtonMethodModule>
("FSHOUnSetup"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<GReKOUpdateSol, NewtonIteratorData, NewtonMethodModule>
("GReKOUpdateSol"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<ImposeHSEquilibriumUpdateSol, NewtonIteratorData, NewtonMethodModule>("ImposeHSEquilibriumUpdateSol"));

Factory<ConvergenceMethod>::getInstance().regist
   (new Environment::ObjectProvider<Linearized, ConvergenceMethod, NewtonMethodModule, 1>
   ("Linearized"));

Factory<ConvergenceMethod>::getInstance().regist
   (new Environment::ObjectProvider<LinearizedBDF2, ConvergenceMethod, NewtonMethodModule, 1>
   ("LinearizedBDF2"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<LinearizedBDF2Setup, NewtonIteratorData, NewtonMethodModule>("LinearizedBDF2Setup"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<LinearizedIntermediate, NewtonIteratorData, NewtonMethodModule>("LinearizedIntermediate"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<LinearizedIntermediateLim, NewtonIteratorData, NewtonMethodModule>("LinearizedIntermediateLim"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<LinearizedPrepare, NewtonIteratorData, NewtonMethodModule>("LinearizedPrepare"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<LinearizedSetup, NewtonIteratorData, NewtonMethodModule>("LinearizedSetup"));

Factory<NewtonIteratorCom>::getInstance().regist
  (new MethodCommandProvider<LinearizedUnSetup, NewtonIteratorData, NewtonMethodModule>("LinearizedUnSetup"));
 
Factory<ConvergenceMethod>::getInstance().regist
   (new Environment::ObjectProvider<NewmarkExplicit, ConvergenceMethod, NewtonMethodModule, 1>
   ("NewmarkExplicit"));
 
Factory<NewtonIteratorCom>::getInstance().regist
  (new MethodCommandProvider<NewmarkExplicitUpdateSol, NewtonIteratorData, NewtonMethodModule> ("NewmarkExplicitUpdateSol"));

 Factory<ConvergenceMethod>::getInstance().regist
   (new Environment::ObjectProvider<NewmarkImplicit, ConvergenceMethod, NewtonMethodModule, 1>
    ("NewmarkImplicit"));

Factory<NewtonIteratorCom>::getInstance().regist
  (new MethodCommandProvider<NewmarkImplicitUpdateSol, NewtonIteratorData, NewtonMethodModule> ("NewmarkImplicitUpdateSol"));

 Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<NewmarkPrepare,NewtonIteratorData, NewtonMethodModule>
    ("NewmarkPrepare"));

Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<NewmarkResetSystem,NewtonIteratorData, NewtonMethodModule>
    ("NewmarkResetSystem"));
 
 Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<NewmarkSetup,NewtonIteratorData, NewtonMethodModule>
    ("NewmarkSetup"));
 
 Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<NewmarkUnSetup,NewtonIteratorData, NewtonMethodModule>
    ("NewmarkUnSetup"));

 Factory<ConvergenceMethod>::getInstance().regist
   (new Environment::ObjectProvider<NewtonIterator, ConvergenceMethod, NewtonMethodModule, 1>
    ("NewtonIterator"));
 
 Factory<ConvergenceMethod>::getInstance().regist
   (new Environment::ObjectProvider<NewtonIteratorCoupling, ConvergenceMethod, NewtonMethodModule, 1>
    ("NewtonIteratorCoupling"));
 
 Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<NullMethodCommand<NewtonIteratorData>,NewtonIteratorData, NewtonMethodModule>("Null"));
 
 Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<ResetSystem, NewtonIteratorData, NewtonMethodModule>
    ("ResetSystem"));

 Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<SelfAdjustUpdateSol, NewtonIteratorData, NewtonMethodModule> 
    ("SelfAdjustUpdateSol"));
 
Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::StdPrepare, NewtonIteratorData, NewtonMethodModule>
("StdPrepare"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::StdSetup, NewtonIteratorData, NewtonMethodModule>
("StdSetup"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::StdUnSetup, NewtonIteratorData, NewtonMethodModule>
("StdUnSetup"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::StdUpdateSol, NewtonIteratorData, NewtonMethodModule>
("StdUpdateSol"));

Factory<NewtonIteratorCom>::getInstance().regist
(new MethodCommandProvider<TurbUpdateSol, NewtonIteratorData, NewtonMethodModule>
("TurbUpdateSol"));

 Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::TwoLayerPrepare, NewtonIteratorData, NewtonMethodModule>("TwoLayerPrepare"));

Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::TwoLayerSetup, NewtonIteratorData, NewtonMethodModule>("TwoLayerSetup"));

Factory<NewtonIteratorCom>::getInstance().regist
  (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::TwoLayerUnSetup, NewtonIteratorData, NewtonMethodModule>("TwoLayerUnSetup"));

Factory<NewtonIteratorCom>::getInstance().regist
  (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::TwoLayerUpdateSol, NewtonIteratorData, NewtonMethodModule>("TwoLayerUpdateSol"));
 
Factory<NewtonIteratorCom>::getInstance().regist
  (new MethodCommandProvider<UpdateSolCoupling, NewtonIteratorData, NewtonMethodModule> 
    ("UpdateSolCoupling"));

 Factory<NewtonIteratorCom>::getInstance().regist
  (new MethodCommandProvider<UpdateSolFVMCC, NewtonIteratorData, NewtonMethodModule> 
    ("UpdateSolFVMCC"));

 Factory<NewtonIteratorCom>::getInstance().regist
   (new MethodCommandProvider<UpdateSolMHD, NewtonIteratorData, NewtonMethodMHDModule> 
    ("UpdateSolMHD"));

 using namespace Petsc;
 
 Factory<ShellPreconditioner>::getInstance().regist
   (new MethodStrategyProvider<BSORPreconditioner, PetscLSSData, ShellPreconditioner,
			       PetscModule>("BSOR"));
 
 Factory<ShellPreconditioner>::getInstance().regist
   (new MethodStrategyProvider<BlockJacobiPreconditioner, PetscLSSData, ShellPreconditioner,
			       PetscModule>("BJacobi"));
 
 Factory<ShellPreconditioner>::getInstance().regist
   (new MethodStrategyProvider<DPLURPreconditioner, PetscLSSData, ShellPreconditioner,
			       PetscModule>("DPLUR"));
 
 Factory<ShellPreconditioner>::getInstance().regist
   (new MethodStrategyProvider<ILUPreconditioner, PetscLSSData, ShellPreconditioner,
			       PetscModule>("ILU"));

 Factory<ShellPreconditioner>::getInstance().regist
   (new MethodStrategyProvider<LUSGSPreconditioner, PetscLSSData, ShellPreconditioner,
			       PetscModule>("LUSGS"));
 
  Factory<PetscLSSCom>::getInstance().regist
  (new MethodCommandProvider<NewParSetup, PetscLSSData, PetscModule>("NewParSetup"));

Factory<ShellPreconditioner>::getInstance().regist
(new MethodStrategyProvider<NullPreconditioner, PetscLSSData, ShellPreconditioner, PetscModule>("Null"));

Factory<PetscLSSCom>::getInstance().regist
    (new MethodCommandProvider<ParJFSetup, PetscLSSData, PetscModule>("ParJFSetup"));
  
  Factory<PetscLSSCom>::getInstance().regist
    (new MethodCommandProvider<ParJFSetupGMRESR, PetscLSSData, PetscModule>("ParJFSetupGMRESR"));
  
  Factory<PetscLSSCom>::getInstance().regist
    (new MethodCommandProvider<ParJFSolveSys, PetscLSSData, PetscModule>("ParJFSolveSys"));
  
  Factory<PetscLSSCom>::getInstance().regist
    (new MethodCommandProvider<ParJFSolveSys, PetscLSSData, PetscModule>("SeqJFSolveSys"));
  
  Factory<PetscLSSCom>::getInstance().regist
    (new MethodCommandProvider<ParJFSolveSysGMRESR, PetscLSSData, PetscModule>("ParJFSolveSysGMRESR"));
  
 Factory<PetscLSSCom>::getInstance().regist
    (new MethodCommandProvider<ParJFSolveSysGMRESR, PetscLSSData, PetscModule>("SeqJFSolveSysGMRESR"));

Factory<PetscLSSCom>::getInstance().regist
    (new MethodCommandProvider<ParMFSetup, PetscLSSData, PetscModule>("ParMFSetup"));
  
Factory<PetscLSSCom>::getInstance().regist
    (new MethodCommandProvider<ParMFSolveSys, PetscLSSData, PetscModule>("ParMFSolveSys"));
  
Factory<PetscLSSCom>::getInstance().regist
    (new MethodCommandProvider<ParMFSolveSys, PetscLSSData, PetscModule>("SeqMFSolveSys"));

Factory<LinearSystemSolver>::getInstance().regist
(new Environment::ObjectProvider<PetscLSS, LinearSystemSolver, PetscModule,1>("PETSC"));

Factory<PetscLSSCom>::getInstance().regist
(new MethodCommandProvider<NullMethodCommand<PetscLSSData>, PetscLSSData, PetscModule>
("Null"));

Factory<PetscLSSCom>::getInstance().regist
(new MethodCommandProvider<StdParSolveSys, PetscLSSData, PetscModule>("StdParSolveSys"));

Factory<PetscLSSCom>::getInstance().regist
(new MethodCommandProvider<StdParSolveSys, PetscLSSData, PetscModule>("StdSeqSolveSys"));

Factory<PetscLSSCom>::getInstance().regist
(new MethodCommandProvider<StdParUnSetup, PetscLSSData, PetscModule>("StdParUnSetup"));

Factory<PetscLSSCom>::getInstance().regist
(new MethodCommandProvider<StdSeqSetup, PetscLSSData, PetscModule>("StdSeqSetup"));

Factory<PetscLSSCom>::getInstance().regist
(new MethodCommandProvider<StdSeqUnSetup, PetscLSSData, PetscModule>("StdSeqUnSetup"));
 
 Factory<ShellPreconditioner>::getInstance().regist
   (new MethodStrategyProvider<TridiagPreconditioner, PetscLSSData, ShellPreconditioner,
			       PetscModule>("Tridiag"));
 
 Factory<PetscLSSCom>::getInstance().regist
   (new MethodCommandProvider<TwoLayerParSetup, PetscLSSData, PetscModule>("TwoLayerParSetup"));
 
 Factory<PetscLSSCom>::getInstance().regist
   (new MethodCommandProvider<TwoLayerParSolveSys, PetscLSSData, PetscModule>
    ("TwoLayerParSolveSys"));
 
 Factory<PetscLSSCom>::getInstance().regist
   (new MethodCommandProvider<TwoLayerSeqSetup, PetscLSSData, PetscModule>
    ("TwoLayerSeqSetup"));

 Factory<PetscLSSCom>::getInstance().regist
   (new MethodCommandProvider<TwoLayerSeqSolveSys, PetscLSSData, PetscModule>
    ("TwoLayerSeqSolveSys"));
}
  
};
  
//////////////////////////////////////////////////////////////////////////////

}
