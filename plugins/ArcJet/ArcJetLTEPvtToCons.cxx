#include "ArcJet/ArcJet.hh"
#include "ArcJetLTEPvtToCons.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "LTE/EulerPvtLTEToCons.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LTE;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////

// reactive gas in LTE 
Environment::ObjectProvider<ArcJetLTEPvtToCons<EulerPvtLTEToCons>, VarSetTransformer, ArcJetModule,1>
arcJetLTE2DPuvtToConsProvider("ArcJetLTE2DPuvtToCons");

Environment::ObjectProvider<ArcJetLTEPvtToCons<EulerPvtLTEToCons>, VarSetTransformer, ArcJetModule,1>
arcJetLTE3DPvtToConsProvider("ArcJetLTE3DPvtToCons");

//////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
