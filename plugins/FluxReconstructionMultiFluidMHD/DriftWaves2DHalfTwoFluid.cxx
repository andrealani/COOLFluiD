#include "FluxReconstructionMultiFluidMHD/DriftWaves2DHalfTwoFluid.hh"
#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMethod/StdSourceTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Maxwell;
using namespace COOLFluiD::Physics::MultiFluidMHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodCommandProvider<DriftWaves2DHalfTwoFluid<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
		       FluxReconstructionSolverData,
		       FluxReconstructionMultiFluidMHDModule>
DriftWaves2DHalfTwoFluidProvider("DriftWaves2DHalfTwoFluid");

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
