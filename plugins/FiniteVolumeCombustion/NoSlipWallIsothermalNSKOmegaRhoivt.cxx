#include "FiniteVolumeCombustion/FiniteVolumeCombustion.hh"
//#include "FiniteVolumeTurb/FiniteVolumeKOmega.hh"
#include "FiniteVolumeCombustion/NoSlipWallIsothermalNSKOmegaRhoivt.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/NSTurbTerm.hh"
#include "NavierStokes/EulerVarSet.hh"     //to be removed?? (used no more)
#include "NavierStokes/Euler2DVarSet.hh"     //added
#include "NavierStokes/Euler3DVarSet.hh"     //added
#include "NavierStokes/MultiScalarVarSet.hh" //added
//#include "NavierStokes/NavierStokesTurbVarSet.hh" //removed (used no more)
#include "NavierStokes/NavierStokes2DVarSet.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"
#include "NEQ/NavierStokesCNEQVarSet.hh"            //added
#include "KOmega/NavierStokesKOmegaVarSet.hh"       //added
#include "NEQKOmega/NavierStokesNEQKOmegaRhoivt.hh" //added



//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NEQKOmega;
using namespace COOLFluiD::Physics::NEQ;
using namespace COOLFluiD::Physics::KOmega;


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NoSlipWallIsothermalNSKOmegaRhoivt
 		      <MultiScalarVarSet
		       <Euler2DVarSet>, 
		       NavierStokesNEQKOmegaRhoivt
		       <NavierStokesKOmegaVarSet<NavierStokesCNEQVarSet<NavierStokes2DVarSet>, 2> > >,
		      CellCenterFVMData,
		      FiniteVolumeCombustionModule>
noSlipWallIsothermalNSKOmegaRhoivt2DFVMCCProvider("NoSlipWallIsothermalNSKOmegaRhoivt2DFVMCC");
      
MethodCommandProvider<NoSlipWallIsothermalNSKOmegaRhoivt
 		      <MultiScalarVarSet<Euler3DVarSet>, 
		       NavierStokesNEQKOmegaRhoivt
		       <NavierStokesKOmegaVarSet<NavierStokesCNEQVarSet<NavierStokes3DVarSet>, 2> > >,
		      CellCenterFVMData,
		      FiniteVolumeCombustionModule>
noSlipWallIsothermalNSKOmegaRhoivt3DFVMCCProvider("NoSlipWallIsothermalNSKOmegaRhoivt3DFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
