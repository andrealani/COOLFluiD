#include "FiniteVolumeGReKO/DistanceBasedExtrapolatorGMoveGReKO.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeGReKO/FiniteVolumeGReKO.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
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

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveGReKO
                       <MultiScalarVarSet<Euler2DVarSet>,
			NavierStokesTurb2DVarSet>,
		       CellCenterFVMData,
		       NodalStatesExtrapolator<CellCenterFVMData>,
		       FiniteVolumeGReKOModule>
distanceBasedExtrapolatorGMoveGReKO2DProvider("DistanceBasedGMoveGReKO2Dnew");
      
//MethodStrategyProvider<DistanceBasedExtrapolatorGMoveGReKO
//                       <MultiScalarVarSet<Euler3DVarSet>,
//                       NavierStokesTurbVarSet>,
//                       CellCenterFVMData,
//                       NodalStatesExtrapolator,
//                       FiniteVolumeGReKOModule>
//distanceBasedExtrapolatorGMoveGReKO3DProvider("DistanceBasedGMoveGReKO3Dnew");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
