#include "HyperPoisson/HyperPoisson.hh"
#include "HyperPoissonPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace HyperPoisson {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<HyperPoissonPhysicalModel<DIM_1D>, PhysicalModelImpl,HyperPoissonModule, 1>
HyperPoisson1DModelProvider("HyperPoisson1D");

Environment::ObjectProvider<HyperPoissonPhysicalModel<DIM_2D>, PhysicalModelImpl,HyperPoissonModule, 1>
HyperPoisson2DModelProvider("HyperPoisson2D");

Environment::ObjectProvider<HyperPoissonPhysicalModel<DIM_3D>, PhysicalModelImpl,HyperPoissonModule, 1>
HyperPoisson3DModelProvider("HyperPoisson3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace HyperPoisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

