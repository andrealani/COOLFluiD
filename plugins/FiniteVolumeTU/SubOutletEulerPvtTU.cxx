#include "FiniteVolumeTU/FiniteVolumeTU.hh"
#include "FiniteVolumeNavierStokes/SubOutletEulerPvt.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler3DRotationVarSet.hh"
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

MethodCommandProvider<SubOutletEulerPvt<Euler3DRotationVarSet>,
                      CellCenterFVMData,
                      FiniteVolumeTUModule>
subOutletEuler3DRotationPvtFVMCCProvider("SubOutletEuler3DRotationPvtFVMCC");
                                          
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics
  
} // namespace COOLFluiD
