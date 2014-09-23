#include "NEQ/NEQ.hh"
#include "NEQ/NavierStokesNEQ.hh"
#include "NEQ/NEQReactionTerm.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokesNEQ<DIM_2D, MultiScalarTerm<NavierStokes::EulerTerm>,
					    NavierStokes::NSTerm, NEQReactionTerm>, PhysicalModelImpl,NEQModule, 1>
navierStokes2DNEQProvider("NavierStokes2DNEQ");
  
Environment::ObjectProvider<NavierStokesNEQ<DIM_3D, MultiScalarTerm<NavierStokes::EulerTerm>,
					    NavierStokes::NSTerm, NEQReactionTerm>, PhysicalModelImpl,NEQModule, 1>
navierStokes3DNEQProvider("NavierStokes3DNEQ");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
