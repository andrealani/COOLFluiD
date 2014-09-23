#include "Framework/DistanceBasedExtrapolator.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DistanceBasedExtrapolator<CellCenterFVMData>,
                       CellCenterFVMData,
                       NodalStatesExtrapolator<CellCenterFVMData>,
                       FiniteVolumeModule>
distanceBasedExtrapolatorProvider("DistanceBased");

//////////////////////////////////////////////////////////////////////////////
  
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
