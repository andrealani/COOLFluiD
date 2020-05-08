#include "KOmega.hh"
#include "NavierStokesKLogOmegaPhysicalModel.hh"
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

Environment::ObjectProvider<NavierStokesKLogOmegaPhysicalModel<DIM_2D>, 
	       PhysicalModelImpl, KOmegaModule, 1>
navierStokes2DKLogOmegaProvider("NavierStokes2DKLogOmega");

Environment::ObjectProvider<NavierStokesKLogOmegaPhysicalModel<DIM_3D>, 
	       PhysicalModelImpl, KOmegaModule, 1>
navierStokes3DKLogOmegaProvider("NavierStokes3DKLogOmega");

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

