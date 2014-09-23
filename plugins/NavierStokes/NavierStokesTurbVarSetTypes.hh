#ifndef COOLFluiD_Physics_NavierStokes_NavierStokesTurbVarSetTypes_hh
#define COOLFluiD_Physics_NavierStokes_NavierStokesTurbVarSetTypes_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesTurbVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////
	
namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

typedef  NavierStokesTurbVarSet<NavierStokes2DVarSet, 0> NavierStokesTurb2DVarSet;
typedef  NavierStokesTurbVarSet<NavierStokes3DVarSet, 0> NavierStokesTurb3DVarSet;

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes
    
  } // namespace Physics
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokesTurbVarSetTypes_hh
