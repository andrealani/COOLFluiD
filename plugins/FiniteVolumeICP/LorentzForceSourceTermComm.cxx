#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/BaseTerm.hh"
#include "FiniteVolumeICP/FiniteVolumeICP.hh"
#include "FiniteVolumeICP/LorentzForceSourceTermComm.hh"
#include "ICP/ICPReactionTerm.hh"

/////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::ICP;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Numerics {

      namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LorentzForceSourceTerm<ICPReactionTerm<BaseTerm> >,
		      DataProcessingData, FiniteVolumeICPModule>
LorentzForceSourceTermCommProvider("LorentzForceSourceTermComm");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

