#include "ArcJet/ArcJet.hh"
#include "ArcJetConvVarSet.hh"
#include "NavierStokes/Euler2DCons.hh"
#include "NavierStokes/Euler3DCons.hh"
#include "LTE/Euler2DPuvtLTE.hh"
#include "LTE/Euler3DPvtLTE.hh"
#include "NavierStokes/Euler2DPuvt.hh"
#include "NavierStokes/Euler3DPvt.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LTE;
using namespace COOLFluiD::Physics::ArcJet;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

/// Cons VarSet does nothing but is needed for giving variables so it 
/// is parametrized with the non-LTE varset.
Environment::ObjectProvider<ArcJetConvVarSet<Euler2DPuvtLTE>, 
			    ConvectiveVarSet, ArcJetModule, 1> 
comparcJetLTE2DPuvtProvider("ArcJetLTE2DPuvt");

Environment::ObjectProvider<ArcJetConvVarSet<Euler3DPvtLTE>, 
			    ConvectiveVarSet, ArcJetModule, 1> 
comparcJetLTE3DPvtProvider("ArcJetLTE3DPvt");
    
Environment::ObjectProvider<ArcJetConvVarSet<Euler2DPuvt>, 
			    ConvectiveVarSet, ArcJetModule, 1> 
comparcJetPG2DPuvtProvider("ArcJetPG2DPuvt");

Environment::ObjectProvider<ArcJetConvVarSet<Euler3DPvt<Euler3DVarSet> >, 
			    ConvectiveVarSet, ArcJetModule, 1> 
comparcJetPG3DPvtProvider("ArcJetPG3DPvt");
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
