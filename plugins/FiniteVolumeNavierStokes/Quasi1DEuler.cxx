#include "Quasi1DEuler.hh"  
#include "NavierStokes/Euler1DVarSet.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics { 
   
    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<Quasi1DEuler
		       <Euler1DVarSet>,
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
quasi1DEulerSTFVMCCProvider("Quasi1DEulerST");

MethodStrategyProvider<Quasi1DEuler
		       <MultiScalarVarSet<Euler1DVarSet> >,
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
quasi1DEulerMSSTVMCCProvider("Quasi1DEulerMSST");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
