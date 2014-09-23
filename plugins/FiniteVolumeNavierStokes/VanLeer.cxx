#include "FiniteVolumeNavierStokes/VanLeer1D.hh"
#include "FiniteVolumeNavierStokes/VanLeer2D.hh"
#include "FiniteVolumeNavierStokes/VanLeer3D.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<VanLeer1D<Euler1DVarSet>,
			CellCenterFVMData,
			FluxSplitter<CellCenterFVMData>,
			FiniteVolumeNavierStokesModule>
vanLeer1DProvider("VanLeer1D");
  
MethodStrategyProvider<VanLeer1D
		       <MultiScalarVarSet<Euler1DVarSet> >,
			CellCenterFVMData,
			FluxSplitter<CellCenterFVMData>,
			FiniteVolumeNavierStokesModule>
vanLeerMS1DProvider("VanLeerMS1D");

MethodStrategyProvider<VanLeer2D<Euler2DVarSet>,
			CellCenterFVMData,
			FluxSplitter<CellCenterFVMData>,
			FiniteVolumeNavierStokesModule>
vanLeer2DProvider("VanLeer2D");
  
MethodStrategyProvider<VanLeer2D
		       <MultiScalarVarSet<Euler2DVarSet> >,
			CellCenterFVMData,
			FluxSplitter<CellCenterFVMData>,
			FiniteVolumeNavierStokesModule>
vanLeerMS2DProvider("VanLeerMS2D");

MethodStrategyProvider<VanLeer3D<Euler3DVarSet>,
			CellCenterFVMData,
			FluxSplitter<CellCenterFVMData>,
			FiniteVolumeNavierStokesModule>
vanLeer3DProvider("VanLeer3D");
  
MethodStrategyProvider<VanLeer3D
		       <MultiScalarVarSet<Euler3DVarSet> >,
			CellCenterFVMData,
			FluxSplitter<CellCenterFVMData>,
			FiniteVolumeNavierStokesModule>
vanLeerMS3DProvider("VanLeerMS3D");
  
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
