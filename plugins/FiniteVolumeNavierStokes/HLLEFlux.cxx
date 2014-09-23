#include "FiniteVolume/HLLEFlux.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<HLLEFlux<Euler1DVarSet>,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeNavierStokesModule>
hlle1DFluxSplitterProvider("HLLE1D");
     
MethodStrategyProvider<HLLEFlux<Euler2DVarSet>,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeNavierStokesModule>
hlle2DFluxSplitterProvider("HLLE2D");
      
MethodStrategyProvider<HLLEFlux<Euler3DVarSet>,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeNavierStokesModule>
hlle3DFluxSplitterProvider("HLLE3D");
      
MethodStrategyProvider<HLLEFlux<MultiScalarVarSet<Euler1DVarSet> >,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeNavierStokesModule>
hlleMS1DFluxSplitterProvider("HLLEMS1D");
     
MethodStrategyProvider<HLLEFlux<MultiScalarVarSet<Euler2DVarSet> >,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeNavierStokesModule>
hlleMS2DFluxSplitterProvider("HLLEMS2D");
      
MethodStrategyProvider<HLLEFlux<MultiScalarVarSet<Euler3DVarSet> >,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeNavierStokesModule>
hlleMS3DFluxSplitterProvider("HLLEMS3D");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
