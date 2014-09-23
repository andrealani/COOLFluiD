#ifndef COOLFluiD_Physics_KOmega_NavierStokesKOmegaVarSetTypes_hh
#define COOLFluiD_Physics_KOmega_NavierStokesKOmegaVarSetTypes_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesTurbVarSetTypes.hh"
#include "KOmega/NavierStokesKOmegaVarSet.hh"
#include "KOmega/NavierStokesKOmegaPvt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

      //NavierStokesKOmegaPvt<NavierStokesKOmegaSSTVarSet<NavierStokes2DVarSet, 0>
      //NavierStokesKOmegaPvt<NavierStokesKOmegaVarSet<NavierStokes2DVarSet, 0>


//////////////////////////////////////////////////////////////////////////////

typedef NavierStokesKOmegaVarSet<NavierStokes::NavierStokes2DVarSet, 0> NavierStokes2DKOmega;
 
typedef NavierStokesKOmegaVarSet<NavierStokes::NavierStokes3DVarSet, 0> NavierStokes3DKOmega;
 
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

#endif // COOLFluiD_Physics_KOmega_NavierStokesKOmegaVarSetTypes_hh
