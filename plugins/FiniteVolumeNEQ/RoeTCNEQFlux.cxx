#include "FiniteVolume/RoeSAFlux.hh"
#include "FiniteVolume/RoeVLPAFlux.hh"
#include "FiniteVolume/RoeEntropyFixFlux.hh"
#include "FiniteVolumeNEQ/RoeTCNEQFlux.hh"
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

Framework::MethodStrategyProvider<RoeTCNEQFlux<RoeSAFlux,MultiScalarVarSet<Euler2DVarSet> >,
				  CellCenterFVMData,
				  Framework::FluxSplitter<CellCenterFVMData>,
				  FiniteVolumeNEQModule>
roeTCNEQ2DSAFluxProvider("RoeTCNEQ2DSA");
      
Framework::MethodStrategyProvider<RoeTCNEQFlux<RoeVLPAFlux,MultiScalarVarSet<Euler2DVarSet> >,
				  CellCenterFVMData,
				  Framework::FluxSplitter<CellCenterFVMData>,
				  FiniteVolumeNEQModule>
roeTCNEQ2DVLPAFluxProvider("RoeTCNEQ2DVLPA");

Framework::MethodStrategyProvider<RoeTCNEQFlux<RoeFlux,MultiScalarVarSet<Euler2DVarSet> >,
				  CellCenterFVMData,
				  Framework::FluxSplitter<CellCenterFVMData>,
				  FiniteVolumeNEQModule>
roeTCNEQ2DFluxProvider("RoeTCNEQ2D");
      

Framework::MethodStrategyProvider<RoeTCNEQFlux<RoeSAFlux,MultiScalarVarSet<Euler3DVarSet> >,
				  CellCenterFVMData,
				  Framework::FluxSplitter<CellCenterFVMData>,
				  FiniteVolumeNEQModule>
roeTCNEQ3DSAFluxProvider("RoeTCNEQ3DSA");
      
Framework::MethodStrategyProvider<RoeTCNEQFlux<RoeVLPAFlux,MultiScalarVarSet<Euler3DVarSet> >,
				  CellCenterFVMData,
				  Framework::FluxSplitter<CellCenterFVMData>,
				  FiniteVolumeNEQModule>
roeTCNEQ3DVLPAFluxProvider("RoeTCNEQ3DVLPA");

Framework::MethodStrategyProvider<RoeTCNEQFlux<RoeFlux,MultiScalarVarSet<Euler3DVarSet> >,
				  CellCenterFVMData,
				  Framework::FluxSplitter<CellCenterFVMData>,
				  FiniteVolumeNEQModule>
roeTCNEQ3DFluxProvider("RoeTCNEQ3D");

Framework::MethodStrategyProvider<RoeTCNEQFlux<RoeEntropyFixFlux,
					       MultiScalarVarSet<Euler1DVarSet> >,
				  CellCenterFVMData,
				  Framework::FluxSplitter<CellCenterFVMData>,
				  FiniteVolumeNEQModule>
roeTCNEQ1DFluxProvider("RoeEntropyFixTCNEQ1D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
