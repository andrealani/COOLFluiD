#include "NavierStokes/NavierStokes.hh"
#include "NavierStokesPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokesPhysicalModel<DIM_1D>, PhysicalModelImpl, NavierStokesModule, 1>
navierStokes1DProvider("NavierStokes1D");

Environment::ObjectProvider<NavierStokesPhysicalModel<DIM_2D>, PhysicalModelImpl, NavierStokesModule, 1>
navierStokes2DProvider("NavierStokes2D");

Environment::ObjectProvider<NavierStokesPhysicalModel<DIM_3D>, PhysicalModelImpl, NavierStokesModule, 1>
navierStokes3DProvider("NavierStokes3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

