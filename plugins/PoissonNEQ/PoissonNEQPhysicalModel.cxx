#include "PoissonNEQ/PoissonNEQ.hh"
#include "PoissonNEQ/PoissonNEQPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace PoissonNEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<PoissonNEQPhysicalModel<DIM_2D>, PhysicalModelImpl, PoissonNEQModule, 1>
poissonNEQ2DProvider("PoissonNEQ2D");

Environment::ObjectProvider<PoissonNEQPhysicalModel<DIM_3D>, PhysicalModelImpl, PoissonNEQModule, 1>
poissonNEQ3DProvider("PoissonNEQ3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace PoissonNEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
