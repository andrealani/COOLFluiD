#include "KOmega.hh"
#include "NavierStokesKOmegaPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokesKOmegaPhysicalModel<DIM_2D>, 
	       PhysicalModelImpl, KOmegaModule, 1>
navierStokes2DKOmegaProvider("NavierStokes2DKOmega");

Environment::ObjectProvider<NavierStokesKOmegaPhysicalModel<DIM_3D>, 
	       PhysicalModelImpl, KOmegaModule, 1>
navierStokes3DKOmegaProvider("NavierStokes3DKOmega");

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

