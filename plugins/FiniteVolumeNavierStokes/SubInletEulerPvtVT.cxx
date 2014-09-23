#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubInletEulerPvtVT.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/EulerVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubInletEulerPvtVT<EulerVarSet>, 
		      CellCenterFVMData, 
		      FiniteVolumeNavierStokesModule>
subInletEulerPvtVTFVMCCProvider("SubInletEulerPvtVTFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
