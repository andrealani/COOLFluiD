#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/NSFlux.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MHD/MHDProjectionDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NSFlux<MHDProjectionDiffVarSet>,
                       CellCenterFVMData,
		       ComputeDiffusiveFlux,
                       FiniteVolumeMHDModule>
navierStokesMHDDiffusiveFluxProvider("NavierStokesMHD");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
