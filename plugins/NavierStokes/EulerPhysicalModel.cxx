#include "NavierStokes/NavierStokes.hh"
#include "EulerPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerPhysicalModel<DIM_1D>, PhysicalModelImpl,NavierStokesModule, 1>
euler1DModelProvider("Euler1D");

Environment::ObjectProvider<EulerPhysicalModel<DIM_2D>, PhysicalModelImpl,NavierStokesModule, 1>
euler2DModelProvider("Euler2D");

Environment::ObjectProvider<EulerPhysicalModel<DIM_3D>, PhysicalModelImpl,NavierStokesModule, 1>
euler3DModelProvider("Euler3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

