#include "NavierStokes/NavierStokes.hh"
#include "Euler3DPvt.hh"
#include "Euler3DVarSet.hh"
#include "Euler3DRotationVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler3DPvt<Euler3DVarSet>, ConvectiveVarSet, NavierStokesModule, 1> 
euler3DPvtProvider("Euler3DPvt");

Environment::ObjectProvider<Euler3DPvt<Euler3DRotationVarSet>, ConvectiveVarSet, NavierStokesModule, 1> 
euler3DRotationPvtProvider("Euler3DRotationPvt");

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
