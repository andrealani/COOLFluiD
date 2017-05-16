#include "Environment/FileHandlerInputConcrete.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/DirectFileAccess.hh"
#include "Environment/Environment.hh"
#include "Environment/FileHandlerOutputConcrete.hh"
#include "Environment/DirectFileWrite.hh"

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

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

class PluginsRegister {
public:

void registerAll() 
{
using namespace Environment;
using namespace Framework;

Factory<FileHandlerInput>::getInstance().regist
(new Environment::ObjectProvider<Environment::FileHandlerInputConcrete< Environment::DirectFileAccess >, 
				 Environment::FileHandlerInput,
				 Environment::EnvironmentModule >("DirectFileAccess"));

Factory<FileHandlerOutput>::getInstance().regist
(new Environment::ObjectProvider<Environment::FileHandlerOutputConcrete< DirectFileWrite >, 
				 Environment::FileHandlerOutput, 
				 Environment::EnvironmentModule >("DirectFileWrite")); 

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
   
      
    }
  };
  
//////////////////////////////////////////////////////////////////////////////

}
