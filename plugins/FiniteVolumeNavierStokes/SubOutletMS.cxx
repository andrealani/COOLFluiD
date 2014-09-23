#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubOutletMS.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "FiniteVolumeNavierStokes/SubOutletEulerPvt.hh"

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

MethodCommandProvider<SubOutletMS<SubOutletEulerPvt<MultiScalarVarSet<Euler2DVarSet> > >, 
		      CellCenterFVMData, 
		      FiniteVolumeNavierStokesModule>
subOutletPvtMS2DFVMCCProvider("SubOutletMS2DpvtFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics
  
} // namespace COOLFluiD

