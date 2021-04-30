#ifndef COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaVarSetTypes_hh
#define COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaVarSetTypes_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesTurbVarSetTypes.hh"
#include "KOmega/NavierStokesKLogOmegaVarSet.hh"
#include "KOmega/NavierStokesKLogOmegaPvt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

      //NavierStokesKOmegaPvt<NavierStokesKOmegaSSTVarSet<NavierStokes2DVarSet, 0>
      //NavierStokesKOmegaPvt<NavierStokesKOmegaVarSet<NavierStokes2DVarSet, 0>


//////////////////////////////////////////////////////////////////////////////

typedef NavierStokesKLogOmegaVarSet<NavierStokes::NavierStokes2DVarSet, 0> NavierStokes2DKLogOmega;
 
typedef NavierStokesKLogOmegaVarSet<NavierStokes::NavierStokes3DVarSet, 0> NavierStokes3DKLogOmega;

typedef NavierStokesKLogOmegaVarSet<NavierStokes::NavierStokesVarSet, 0> NavierStokesKLogOmega;

 
      /* typedef NavierStokesKOmegaPvt<NavierStokesKOmegaVarSet
	 <NavierStokes::NavierStokes2DVarSet, 0> > NavierStokes2DKOmegaPuvt;
	 
	 typedef NavierStokesKOmegaPvt<NavierStokesKOmegaVarSet
	 <NavierStokes::NavierStokes3DVarSet, 0> > NavierStokes3DKOmegaPvt;
      */  
//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaVarSetTypes_hh
