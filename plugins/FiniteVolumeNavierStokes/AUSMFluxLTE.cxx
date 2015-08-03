#include "Framework/MethodStrategyProvider.hh"

#include "FiniteVolumeNavierStokes/AUSMPlusUpFlux.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeLTE.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "LTE/Euler3DPvtLTE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LTE;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<AUSMPlusUpFlux <MultiScalarVarSet<Euler3DPvtLTE> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeLTEModule>
ausmPlusUpTurb3dLTEProvider("AUSMPlusUpTurb3DLTE");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
