#include "ArcJet/ArcJetSALTE.hh"
#include "ArcJet/ArcJetDiffVarSet.hh"
#include "SA/NavierStokesSAPvtLTE.hh"
#include "SA/NavierStokesSAVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::SA;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////
      
Environment::ObjectProvider<ArcJetDiffVarSet<NavierStokesSAPvtLTE<NavierStokesSAVarSet<NavierStokes2DVarSet, 0> > >,
			    DiffusiveVarSet, 
			    ArcJetSALTEModule, 2> 
arcJetSALTE2DPuvtDiffProvider("ArcJetSALTE2DPuvt");
      
Environment::ObjectProvider<ArcJetDiffVarSet<NavierStokesSAPvtLTE<NavierStokesSAVarSet<NavierStokes3DVarSet, 0> > >,
			    DiffusiveVarSet, 
			    ArcJetSALTEModule, 2> 
arcJetSALTE3DPvtDiffProvider("ArcJetSALTE3DPvt");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
