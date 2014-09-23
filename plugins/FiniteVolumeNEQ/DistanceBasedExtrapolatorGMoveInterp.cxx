#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolumeNEQ/DistanceBasedExtrapolatorGMoveRhoivt.hh"
#include "FiniteVolumeNEQ/DistanceBasedExtrapolatorGMovePivt.hh"
#include "FiniteVolumeNEQ/DistanceBasedExtrapolatorGMoveCat.hh"
#include "FiniteVolumeNEQ/DistanceBasedExtrapolatorGMoveRhoivtLTE.hh"
#include "FiniteVolumeNEQ/DistanceBasedExtrapolatorGMoveInterp.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveInterp
		       <DistanceBasedExtrapolatorGMoveRhoivt>,
		       CellCenterFVMData,
		       NodalStatesExtrapolator<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
distanceBasedExtrapolatorGMoveRhoivtInterpProvider("DistanceBasedGMoveRhoivtInterp");

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveInterp
		       <DistanceBasedExtrapolatorGMovePivt>,
		       CellCenterFVMData,
		       NodalStatesExtrapolator<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
distanceBasedExtrapolatorGMovePivtInterpProvider("DistanceBasedGMovePivtInterp");

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveInterp
		       <DistanceBasedExtrapolatorGMoveCat>,
		       CellCenterFVMData,
		       NodalStatesExtrapolator<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
distanceBasedExtrapolatorGMoveCatInterpProvider("DistanceBasedGMoveCatInterp");

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveInterp
                       <DistanceBasedExtrapolatorGMoveRhoivtLTE>,
                       CellCenterFVMData,
                       NodalStatesExtrapolator<CellCenterFVMData>,
                       FiniteVolumeNEQModule>
distanceBasedExtrapolatorGMoveRhoivtLTEInterpProvider("DistanceBasedGMoveRhoivtLTEInterp");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
