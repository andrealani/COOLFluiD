#include "LTE.hh"
#include "EulerLTEDemix.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerLTEDemix<DIM_2D, NavierStokes::EulerTerm>, PhysicalModelImpl,LTEModule, 1>
euler2DLTEDemixProvider("Euler2DLTEDemix");

Environment::ObjectProvider<EulerLTEDemix<DIM_3D, NavierStokes::EulerTerm>, PhysicalModelImpl,LTEModule, 1>
euler3DLTEDemixProvider("Euler3DLTEDemix");

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

