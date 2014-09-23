#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubOutletEulerPvt.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
	
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

MethodCommandProvider<SubOutletEulerPvt<Euler2DVarSet>, 
		      CellCenterFVMData, 
		      FiniteVolumeNavierStokesModule>
subOutletEuler2DPvtFVMCCProvider("SubOutletEuler2DPvtFVMCC");

MethodCommandProvider<SubOutletEulerPvt<Euler2DVarSet>, 
		      CellCenterFVMData, 
		      FiniteVolumeNavierStokesModule>
subOutletEuler2DPuvtFVMCCProvider("SubOutletEuler2DPuvtFVMCC");

MethodCommandProvider<SubOutletEulerPvt<Euler3DVarSet>,
                      CellCenterFVMData,
                      FiniteVolumeNavierStokesModule>
subOutletEuler3DPvtFVMCCProvider("SubOutletEuler3DPvtFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics
  
} // namespace COOLFluiD

