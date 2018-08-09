#include "Framework/MethodStrategyProvider.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"

#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/AUSMFluxMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/AUSMPlusUpFluxMultiFluidMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::Maxwell;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace std;
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// AUSM+up

Framework::MethodStrategyProvider<AUSMPlusUpFluxMultiFluid<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
		       FluxReconstructionSolverData,
		       RiemannFlux,
		       FluxReconstructionMultiFluidMHDModule> 
ausmPlusUp2DMultiFluidFRProvider("AUSMPlusUpFluxMultiFluidMHD2D");

Framework::MethodStrategyProvider<AUSMPlusUpFluxMultiFluid<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
		       FluxReconstructionSolverData,
		       RiemannFlux,
		       FluxReconstructionMultiFluidMHDModule>
ausmPlusUp3DMultiFluidFRProvider("AUSMPlusUpFluxMultiFluidMHD3D");

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
