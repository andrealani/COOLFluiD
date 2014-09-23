#include "FiniteVolumeICP/FiniteVolumeICPNEQ.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolumeICP/FiniteVolumeICP.hh"
#include "FiniteVolumeICP/RMSJouleHeatSourceCoupling.hh"
#include "ICP/ICPReactionTerm.hh"
#include "NEQ/NEQReactionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::ICP;
using namespace COOLFluiD::Physics::NEQ;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RMSJouleHeatSourceCoupling<ICPReactionTerm<NEQReactionTerm> >, 
		      DataProcessingData, FiniteVolumeICPNEQModule>
rmsJouleHeatSourceCouplingNEQProvider("RMSJouleHeatSourceCouplingNEQ");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
