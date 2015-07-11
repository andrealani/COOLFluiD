#include "ArcJet/ArcJet.hh"
#include "ArcJet/ArcJetInductionConvVarSet.hh"
#include "NavierStokes/Euler2DCons.hh"
#include "NavierStokes/Euler3DCons.hh"
#include "NavierStokes/Euler2DPuvt.hh"
#include "NavierStokes/Euler3DPvt.hh"
#include "LTE/Euler2DPuvtLTE.hh"
#include "LTE/Euler3DPvtLTE.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LTE;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ArcJetInductionConvVarSet<Euler2DCons>, 
			    ConvectiveVarSet, ArcJetModule, 1> 
arcJetLTE2DConsConvProvider("ArcJet2DConsLTE");
    
Environment::ObjectProvider<ArcJetInductionConvVarSet<Euler2DPuvtLTE>, 
			       ConvectiveVarSet, ArcJetModule, 1> 
arcJetLTE2DPuvtConvProvider("ArcJet2DPuvtLTE");

Environment::ObjectProvider<ArcJetInductionConvVarSet<Euler3DCons>, 
			    ConvectiveVarSet, ArcJetModule, 1> 
arcJetLTE3DConsConvProvider("ArcJet3DConsLTE");
    
Environment::ObjectProvider<ArcJetInductionConvVarSet<Euler3DPvtLTE>, 
			       ConvectiveVarSet, ArcJetModule, 1> 
arcJetLTE3DPvtConvProvider("ArcJet3DPvtLTE");

Environment::ObjectProvider<ArcJetInductionConvVarSet<Euler2DCons>, 
			    ConvectiveVarSet, ArcJetModule, 1> 
arcJet2DConsConvProvider("ArcJet2DCons");
    
Environment::ObjectProvider<ArcJetInductionConvVarSet<Euler2DPuvt>, 
			       ConvectiveVarSet, ArcJetModule, 1> 
arcJet2DPuvtConvProvider("ArcJet2DPuvt");

Environment::ObjectProvider<ArcJetInductionConvVarSet<Euler3DCons>, 
			    ConvectiveVarSet, ArcJetModule, 1> 
arcJet3DConsConvProvider("ArcJet3DCons");
    
Environment::ObjectProvider<ArcJetInductionConvVarSet<Euler3DPvt<Euler3DVarSet> >, 
			    ConvectiveVarSet, ArcJetModule, 1> 
arcJet3DPvtConvProvider("ArcJet3DPvt");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
