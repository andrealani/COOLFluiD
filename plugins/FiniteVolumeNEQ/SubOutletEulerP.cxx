#include "FiniteVolumeNEQ/SubOutletEulerP.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubOutletEulerP<Euler1DVarSet>,
                      CellCenterFVMData,
		      FiniteVolumeNEQModule>
subOutletEuler1DPFVMCCProvider("SubOutletEuler1DPFVMCC");

MethodCommandProvider<SubOutletEulerP<Euler2DVarSet>,
                      CellCenterFVMData,
		      FiniteVolumeNEQModule>
subOutletEuler2DPFVMCCProvider("SubOutletEuler2DPFVMCC");

MethodCommandProvider<SubOutletEulerP<Euler3DVarSet>,
                      CellCenterFVMData,
		      FiniteVolumeNEQModule>
subOutletEuler3DPFVMCCProvider("SubOutletEuler3DPFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
