#include "FiniteVolume/RoeSAFlux.hh"
#include "FiniteVolume/RoeVLPAFlux.hh"
#include "FiniteVolume/RoeEntropyFixFlux.hh"
#include "FiniteVolumeNEQ/RoeVinokurTCNEQFlux.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<RoeVinokurTCNEQFlux<RoeSAFlux,MultiScalarVarSet<Euler2DVarSet> >,
				  CellCenterFVMData,
				  Framework::FluxSplitter<CellCenterFVMData>,
				  FiniteVolumeNEQModule>
roeVinokurTCNEQ2DSAFluxProvider("RoeVinokurTCNEQ2DSA");
      
Framework::MethodStrategyProvider<RoeVinokurTCNEQFlux<RoeVLPAFlux,MultiScalarVarSet<Euler2DVarSet> >,
				  CellCenterFVMData,
				  Framework::FluxSplitter<CellCenterFVMData>,
				  FiniteVolumeNEQModule>
roeVinokurTCNEQ2DVLPAFluxProvider("RoeVinokurTCNEQ2DVLPA");

Framework::MethodStrategyProvider<RoeVinokurTCNEQFlux<RoeFlux,MultiScalarVarSet<Euler2DVarSet> >,
				  CellCenterFVMData,
				  Framework::FluxSplitter<CellCenterFVMData>,
				  FiniteVolumeNEQModule>
roeVinokurTCNEQ2DFluxProvider("RoeVinokurTCNEQ2D");
      

// Framework::MethodStrategyProvider<RoeVinokurTCNEQFlux<RoeSAFlux,MultiScalarVarSet<Euler3DVarSet> >,
// 				  CellCenterFVMData,
// 				  Framework::FluxSplitter<CellCenterFVMData>,
// 				  FiniteVolumeNEQModule>
// roeVinokurTCNEQ3DSAFluxProvider("RoeVinokurTCNEQ3DSA");
// 
// Framework::MethodStrategyProvider<RoeVinokurTCNEQFlux<RoeVLPAFlux,MultiScalarVarSet<Euler3DVarSet> >,
// 				  CellCenterFVMData,
// 				  Framework::FluxSplitter<CellCenterFVMData>,
// 				  FiniteVolumeNEQModule>
// roeVinokurTCNEQ3DVLPAFluxProvider("RoeVinokurTCNEQ3DVLPA");
// 
// Framework::MethodStrategyProvider<RoeVinokurTCNEQFlux<RoeFlux,MultiScalarVarSet<Euler3DVarSet> >,
// 				  CellCenterFVMData,
// 				  Framework::FluxSplitter<CellCenterFVMData>,
// 				  FiniteVolumeNEQModule>
// roeVinokurTCNEQ3DFluxProvider("RoeVinokurTCNEQ3D");

Framework::MethodStrategyProvider<RoeVinokurTCNEQFlux<RoeEntropyFixFlux,
					       MultiScalarVarSet<Euler1DVarSet> >,
				  CellCenterFVMData,
				  Framework::FluxSplitter<CellCenterFVMData>,
				  FiniteVolumeNEQModule>
roeVinokurTCNEQ1DFluxProvider("RoeVinokurEntropyFixTCNEQ1D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
