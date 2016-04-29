#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "FiniteVolume/MeshFittingAlgorithm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace COOLFluiD::Physics::Maxwell;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MeshFittingAlgorithm<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >, 
		      DataProcessingData, 
		      FiniteVolumeMultiFluidMHDModule>
meshFittingAlgorithmMultiFluidMHD2DProvider("MeshFittingAlgorithmMFMHD2D");
      
MethodCommandProvider<MeshFittingAlgorithm<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >, 
		      DataProcessingData, 
		      FiniteVolumeMultiFluidMHDModule>
meshFittingAlgorithmMultiFluidMHD3DProvider("MeshFittingAlgorithmMFMHD3D");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
