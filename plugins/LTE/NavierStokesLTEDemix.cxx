#include "LTE.hh"
#include "NavierStokesLTEDemix.hh"
#include "Environment/ObjectProvider.hh" 
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokesLTEDemix<DIM_2D, EulerTerm>, 
	PhysicalModelImpl, LTEModule, 1> navierStokes2DLTEDemixProvider("NavierStokes2DLTEDemix");

Environment::ObjectProvider<NavierStokesLTEDemix<DIM_3D, EulerTerm>, 
	PhysicalModelImpl, LTEModule, 1> navierStokes3DLTEDemixProvider("NavierStokes3DLTEDemix");

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

