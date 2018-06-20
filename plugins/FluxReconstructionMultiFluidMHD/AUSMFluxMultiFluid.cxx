#include "Framework/MethodStrategyProvider.hh"

#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::Maxwell;
using namespace COOLFluiD::Physics::MultiFluidMHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// AUSM+up

/*MethodStrategyProvider<AUSMPlusUpFluxMultiFluid<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
		       FluxReconstructionSolverData,
		       RiemannFlux,
		       FluxReconstructionNavierStokesModule > 
		       ausmPlusUp2DMultiFluidFRProvider("AUSMPlusUpMultiFluid2D");*/

/*MethodStrategyProvider<AUSMPlusUpFluxMultiFluid<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
		       FluxReconstructionSolverData,
		       RiemannFlux,
		       FluxReconstructionNavierStokesModule > 
		       ausmPlusUp3DMultiFluidFRProvider("AUSMPlusUpMultiFluid3D");*/

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
