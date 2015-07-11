#include "ArcJet/ArcJet.hh"
#include "ArcJetPvtToCons.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/Euler2DPuvtToCons.hh"
#include "NavierStokes/Euler3DPvtToCons.hh"
#include "LTE/Euler2DPuvtLTEToCons.hh"
#include "LTE/Euler3DPvtLTEToCons.hh"

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
Environment::ObjectProvider<ArcJetPvtToCons<Euler2DPuvtToCons>, VarSetTransformer, ArcJetModule,1>
arcJet2DPuvtToConsProvider("ArcJet2DPuvtToCons");

Environment::ObjectProvider<ArcJetPvtToCons<Euler3DPvtToCons>, VarSetTransformer, ArcJetModule,1>
arcJet3DPvtToConsProvider("ArcJet3DPvtToCons");

// reactive gas in LTE 
Environment::ObjectProvider<ArcJetPvtToCons<Euler2DPuvtLTEToCons>, VarSetTransformer, ArcJetModule,1>
arcJet2DPuvtLTEToConsProvider("ArcJet2DPuvtLTEToCons");

Environment::ObjectProvider<ArcJetPvtToCons<Euler3DPvtLTEToCons>, VarSetTransformer, ArcJetModule,1>
arcJet3DPvtLTEToConsProvider("ArcJet3DPvtLTEToCons");

//////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
