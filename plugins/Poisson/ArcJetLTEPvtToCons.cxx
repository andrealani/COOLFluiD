#include "ArcJet/ArcJet.hh"
#include "ArcJetLTEPvtToCons.hh"
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

// reactive gas in LTE 
Environment::ObjectProvider<ArcJetLTEPvtToCons<Euler2DPuvtLTEToCons>, VarSetTransformer, ArcJetModule,1>
arcJetLTE2DPuvtToConsProvider("ArcJetLTE2DPuvtToCons");

Environment::ObjectProvider<ArcJetLTEPvtToCons<Euler3DPvtLTEToCons>, VarSetTransformer, ArcJetModule,1>
arcJetLTE3DPvtToConsProvider("ArcJetLTE3DPvtToCons");

//////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
