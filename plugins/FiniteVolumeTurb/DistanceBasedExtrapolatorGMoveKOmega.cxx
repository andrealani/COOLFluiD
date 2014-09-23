#include "FiniteVolumeTurb/DistanceBasedExtrapolatorGMoveKOmega.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeTurb/FiniteVolumeKOmega.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSetTypes.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveKOmega
                       <EulerVarSet, NavierStokesTurb2DVarSet>,
                       CellCenterFVMData,
                       NodalStatesExtrapolator<CellCenterFVMData>,
                       FiniteVolumeKOmegaModule>
distanceBasedExtrapolatorGMoveKOmega2DProvider("DistanceBasedGMoveKOmega2Dnew");

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveKOmega
                       <EulerVarSet, NavierStokesTurb3DVarSet>,
                       CellCenterFVMData,
                       NodalStatesExtrapolator<CellCenterFVMData>,
                       FiniteVolumeKOmegaModule>
distanceBasedExtrapolatorGMoveKOmega3DProvider("DistanceBasedGMoveKOmega3Dnew");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
