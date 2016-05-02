#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/NSFluxCoupling.hh"
#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "MultiFluidMHD/DiffMFMHDVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MultiFluidMHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NSFlux<DiffMFMHDVarSet>,
                       CellCenterFVMData,
		       ComputeDiffusiveFlux,
                       FiniteVolumeMultiFluidMHDModule>
navierStokesMFDiffusiveFluxProvider("NavierStokesMF");

MethodStrategyProvider<NSFluxCoupling<DiffMFMHDVarSet>,
                       CellCenterFVMData,
		       ComputeDiffusiveFlux,
                       FiniteVolumeMultiFluidMHDModule>
navierStokesMFCouplingDiffusiveFluxProvider("NavierStokesMFCoupling");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
