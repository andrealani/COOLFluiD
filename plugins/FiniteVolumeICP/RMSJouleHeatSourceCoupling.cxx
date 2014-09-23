#include "FiniteVolumeICP/FiniteVolumeICP.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolumeICP/FiniteVolumeICP.hh"
#include "FiniteVolumeICP/RMSJouleHeatSourceCoupling.hh"
#include "ICP/ICPReactionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::ICP;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RMSJouleHeatSourceCoupling<ICPReactionTerm<BaseTerm> >, 
		      DataProcessingData, FiniteVolumeICPModule>
rmsJouleHeatSourceCouplingProvider("RMSJouleHeatSourceCoupling");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
