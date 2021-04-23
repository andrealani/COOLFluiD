#include "GammaAlpha.hh"
#include "NavierStokesGammaAlphaPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokesGammaAlphaPhysicalModel<DIM_2D>, 
	       PhysicalModelImpl, GammaAlphaModule, 1>
navierStokes2DGammaAlphaProvider("NavierStokes2DGammaAlpha");

Environment::ObjectProvider<NavierStokesGammaAlphaPhysicalModel<DIM_3D>, 
	       PhysicalModelImpl, GammaAlphaModule, 1>
navierStokes3DGammaAlphaProvider("NavierStokes3DGammaAlpha");

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace GammaAlpha

  } // end of namespace Physics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

