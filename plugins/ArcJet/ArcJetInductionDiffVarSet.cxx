#include "ArcJet/ArcJet.hh"
#include "ArcJet/ArcJetInductionDiffVarSet.hh"
#include "ArcJet/ArcJetTerm.hh"
#include "ArcJet/ArcJetInductionTerm.hh"
#include "LTE/NavierStokes2DPuvtLTE.hh"
#include "LTE/NavierStokes3DPvtLTE.hh"
#include "NavierStokes/NavierStokes2DPuvt.hh"
#include "NavierStokes/NavierStokes3DPvt.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/MultiScalarTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LTE;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////
      
Environment::ObjectProvider<ArcJetInductionDiffVarSet
			    <NavierStokes2DPuvtLTE<ArcJetInductionTerm<EulerTerm> >, 
			     ArcJetTerm<Framework::BaseTerm> >, 
			    DiffusiveVarSet, 
			    ArcJetModule, 2> 
arcJet2DPuvtLTEDiffProvider("ArcJet2DPuvtLTE");

Environment::ObjectProvider<ArcJetInductionDiffVarSet
			    <NavierStokes3DPvtLTE, 
			    ArcJetTerm<Framework::BaseTerm> >,
                            DiffusiveVarSet, 
                            ArcJetModule, 2> 
arcJet3DPvtLTEDiffProvider("ArcJet3DPvtLTE");

Environment::ObjectProvider<ArcJetInductionDiffVarSet
			    <NavierStokes2DPuvt, 
			     ArcJetTerm<Framework::BaseTerm> >, 
			    DiffusiveVarSet, 
			    ArcJetModule, 2> 
arcJet2DPuvtDiffProvider("ArcJet2DPuvt");

Environment::ObjectProvider<ArcJetInductionDiffVarSet
			    <NavierStokes3DPvt, 
			    ArcJetTerm<Framework::BaseTerm> >,
                            DiffusiveVarSet, 
                            ArcJetModule, 2> 
arcJet3DPvtDiffProvider("ArcJet3DPvt");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
