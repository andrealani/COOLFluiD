#include "ArcJet/ArcJetSALTE.hh"
#include "ArcJet/ArcJetLTEPvtToCons.hh"
#include "Environment/ObjectProvider.hh"
#include "LTE/EulerPvtLTEToCons.hh"
#include "SA/EulerSAPvtToCons.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LTE;
using namespace COOLFluiD::Physics::SA;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////

// reactive gas in LTE with SA
Environment::ObjectProvider<ArcJetLTEPvtToCons<EulerSAPvtToCons<EulerPvtLTEToCons> >, 
			    VarSetTransformer, ArcJetSALTEModule,1>
arcJetSALTE2DPuvtToConsProvider("ArcJetSALTE2DPuvtToCons");

Environment::ObjectProvider<ArcJetLTEPvtToCons<EulerSAPvtToCons<EulerPvtLTEToCons> >, 
			    VarSetTransformer, ArcJetSALTEModule,1>
arcJetSALTE3DPvtToConsProvider("ArcJetSALTE3DPvtToCons");

//////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
