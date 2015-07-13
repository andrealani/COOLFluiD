#include "ArcJet/ArcJet.hh"
#include "ArcJet/ArcJetPvtToCons.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/EulerPvtToCons.hh"
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

// perfect non reactive gas case
Environment::ObjectProvider<ArcJetPvtToCons<EulerPvtToCons>, VarSetTransformer, ArcJetModule,1>
arcJet2DPuvtToConsProvider("ArcJet2DPuvtToCons");

Environment::ObjectProvider<ArcJetPvtToCons<EulerPvtToCons>, VarSetTransformer, ArcJetModule,1>
arcJet3DPvtToConsProvider("ArcJet3DPvtToCons");

// reactive gas in LTE 
Environment::ObjectProvider<ArcJetPvtToCons<EulerPvtLTEToCons>, VarSetTransformer, ArcJetModule,1>
arcJet2DPuvtLTEToConsProvider("ArcJet2DPuvtLTEToCons");

Environment::ObjectProvider<ArcJetPvtToCons<EulerPvtLTEToCons>, VarSetTransformer, ArcJetModule,1>
arcJet3DPvtLTEToConsProvider("ArcJet3DPvtLTEToCons");

//////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
