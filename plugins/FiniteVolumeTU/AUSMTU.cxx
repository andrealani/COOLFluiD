#include "FiniteVolumeNavierStokes/AUSMPlusFlux.hh"
#include "FiniteVolumeNavierStokes/AUSMPlusUpFlux.hh"
#include "FiniteVolumeNavierStokes/AUSMLiouSteffenFlux.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler3DRotationVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeTU/FiniteVolumeTU.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/// AUSM+up
MethodStrategyProvider<AUSMPlusUpFlux<Euler3DRotationVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeTUModule>
ausmPlusUp3dRotationProvider("PlusUp3DRotation");

MethodStrategyProvider<AUSMPlusUpFlux<MultiScalarVarSet<Euler3DRotationVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeTUModule>
ausmPlusUpMS3dRotationProvider("PlusUpMS3DRotation");

/// AUSM+
MethodStrategyProvider<AUSMPlusFlux<Euler3DRotationVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeTUModule>
ausmPlus3dRotationProvider("Plus3DRotation");

MethodStrategyProvider<AUSMPlusFlux<MultiScalarVarSet<Euler3DRotationVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeTUModule>
ausmPlusMS3dRotationProvider("PlusMS3DRotation");

/// AUSM Liou-Steffen
MethodStrategyProvider<AUSMLiouSteffenFlux<Euler3DRotationVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeTUModule>
ausmLiouSteffen3dRotationProvider("LiouSteffen3DRotation");

MethodStrategyProvider<AUSMLiouSteffenFlux<MultiScalarVarSet<Euler3DRotationVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeTUModule>
ausmLiouSteffenMS3dRotationProvider("LiouSteffenMS3DRotation");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
