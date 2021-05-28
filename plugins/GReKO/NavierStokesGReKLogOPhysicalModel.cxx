#include "GReKO.hh"
#include "NavierStokesGReKLogOPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokesGReKLogOPhysicalModel<DIM_2D>, 
	       PhysicalModelImpl, GReKOModule, 1>
navierStokes2DGReKLogOProvider("NavierStokes2DGReKLogO");

Environment::ObjectProvider<NavierStokesGReKLogOPhysicalModel<DIM_3D>, 
	       PhysicalModelImpl, GReKOModule, 1>
navierStokes3DGReKLogOProvider("NavierStokes3DGReKLogO");

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace GReKO

  } // end of namespace Physics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

