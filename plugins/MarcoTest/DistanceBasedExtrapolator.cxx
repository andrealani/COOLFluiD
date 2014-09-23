#include "Framework/DistanceBasedExtrapolator.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "MarcoTest/MarcoTest.hh"
#include "MarcoTest/MarcoTestMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MarcoTest {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DistanceBasedExtrapolator<MarcoTestMethodData>,
                       MarcoTestMethodData,
                       NodalStatesExtrapolator<MarcoTestMethodData>,
                       MarcoTestModule>
distanceBasedExtrapolatorMarcoTestProvider("DistanceBased");

//////////////////////////////////////////////////////////////////////////////

  } // namespace MarcoTest

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
