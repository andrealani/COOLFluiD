#include "NavierStokesChannelSourceTerm.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NavierStokesChannelSourceTerm<EulerVarSet,
                                                     NavierStokesVarSet>,
		                   CellCenterFVMData, 
		                   ComputeSourceTerm<CellCenterFVMData>,
		                   FiniteVolumeNavierStokesModule>
navierStokes3DChannelSTFVMCCProvider("NavierStokesChannelST");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
