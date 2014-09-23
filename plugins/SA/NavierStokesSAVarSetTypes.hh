#ifndef COOLFluiD_Physics_SA_NavierStokesSAVarSetTypes_hh
#define COOLFluiD_Physics_SA_NavierStokesSAVarSetTypes_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesTurbVarSetTypes.hh"
#include "SA/NavierStokesSAVarSet.hh"
#include "SA/NavierStokesSAPvt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

typedef NavierStokesSAVarSet<NavierStokes::NavierStokes2DVarSet, 0> NavierStokes2DSA;
 
typedef NavierStokesSAVarSet<NavierStokes::NavierStokes3DVarSet, 0> NavierStokes3DSA;

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_SA_NavierStokesSAVarSetTypes_hh
