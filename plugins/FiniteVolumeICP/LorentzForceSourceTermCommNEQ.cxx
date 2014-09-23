#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataProcessing.hh"
#include "FiniteVolumeICP/FiniteVolumeICP.hh"
#include "FiniteVolumeICP/LorentzForceSourceTermComm.hh"
#include "ICP/ICPReactionTerm.hh"
#include "NEQ/NEQReactionTerm.hh"

/////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::ICP;
using namespace COOLFluiD::Physics::NEQ;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Numerics {

      namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LorentzForceSourceTerm<ICPReactionTerm<NEQReactionTerm> >,
		      DataProcessingData, FiniteVolumeICPModule>
LorentzForceSourceTermCommNEQProvider("LorentzForceSourceTermCommNEQ");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

