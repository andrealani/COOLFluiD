#include "GReKO.hh"
#include "NavierStokesGReKOPhysicalModel.hh"
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

Environment::ObjectProvider<NavierStokesGReKOPhysicalModel<DIM_2D>, 
	       PhysicalModelImpl, GReKOModule, 1>
navierStokes2DGReKOProvider("NavierStokes2DGReKO");

//Environment::ObjectProvider<NavierStokesGReKOPhysicalModel<DIM_3D>, 
//	       PhysicalModelImpl, GReKOModule, 1>
//navierStokes3DGReKOProvider("NavierStokes3DGReKO");

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace GReKO

  } // end of namespace Physics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

