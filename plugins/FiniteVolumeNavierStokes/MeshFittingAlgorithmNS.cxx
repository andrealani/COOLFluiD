#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "FiniteVolume/MeshFittingAlgorithm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MeshFittingAlgorithm<EulerVarSet>, 
		      DataProcessingData, 
		      FiniteVolumeNavierStokesModule>
meshFittingAlgorithmNSProvider("MeshFittingAlgorithmNS");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
