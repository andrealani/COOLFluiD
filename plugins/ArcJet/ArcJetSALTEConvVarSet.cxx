#include "ArcJet/ArcJetSALTE.hh"
#include "ArcJet/ArcJetConvVarSet.hh"
#include "LTE/Euler2DPuvtLTE.hh"
#include "LTE/Euler3DPvtLTE.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "SA/EulerSAVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LTE;
using namespace COOLFluiD::Physics::ArcJet;
using namespace COOLFluiD::Physics::SA;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

/// Cons VarSet does nothing but is needed for giving variables so it 
/// is parametrized with the non-LTE varset.
Environment::ObjectProvider<ArcJetConvVarSet<EulerSAVarSet<MultiScalarVarSet<Euler2DPuvtLTE>,0> >, 
			    ConvectiveVarSet, ArcJetSALTEModule, 1> 
arcJetSALTE2DPuvtProvider("ArcJetSALTE2DPuvt");

Environment::ObjectProvider<ArcJetConvVarSet<EulerSAVarSet<MultiScalarVarSet<Euler3DPvtLTE>,0> >, 
			    ConvectiveVarSet, ArcJetSALTEModule, 1> 
arcJetSALTE3DPvtProvider("ArcJetSALTE3DPvt");
        
//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
