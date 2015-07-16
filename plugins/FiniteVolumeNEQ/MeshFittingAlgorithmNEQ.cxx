#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolumeNEQ/ThermNEQST.hh"
#include "FiniteVolumeNEQ/ChemNEQST.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolume/MeshFittingAlgorithm.hh"

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

MethodCommandProvider<MeshFittingAlgorithm<MultiScalarVarSet<Euler1DVarSet> >,
		      DataProcessingData, 
		      FiniteVolumeNEQModule>
meshFittingAlgorithmNEQ1DProvider("MeshFittingAlgorithmNEQ1D");

MethodCommandProvider<MeshFittingAlgorithm<MultiScalarVarSet<Euler2DVarSet> >,
		      DataProcessingData, 
		      FiniteVolumeNEQModule>
meshFittingAlgorithmNEQ2DProvider("MeshFittingAlgorithmNEQ2D");

MethodCommandProvider<MeshFittingAlgorithm<MultiScalarVarSet<Euler3DVarSet> >,
		      DataProcessingData, 
		      FiniteVolumeNEQModule>
meshFittingAlgorithmNEQ3DProvider("MeshFittingAlgorithmNEQ3D");
 
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
