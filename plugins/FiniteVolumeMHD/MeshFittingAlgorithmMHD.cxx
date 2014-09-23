#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MHD/MHD2DVarSet.hh"
#include "MHD/MHD3DVarSet.hh"
#include "MHD/MHD2DProjectionVarSet.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "FiniteVolume/MeshFittingAlgorithm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MeshFittingAlgorithm<MHD2DVarSet>, 
		      DataProcessingData, 
		      FiniteVolumeMHDModule>
meshFittingAlgorithmMHD2DProvider("MeshFittingAlgorithmMHD2D");
      
MethodCommandProvider<MeshFittingAlgorithm<MHD3DVarSet>, 
		      DataProcessingData, 
		      FiniteVolumeMHDModule>
meshFittingAlgorithmMHD3DProvider("MeshFittingAlgorithmMHD3D");

MethodCommandProvider<MeshFittingAlgorithm<MHD2DProjectionVarSet>, 
		      DataProcessingData, 
		      FiniteVolumeMHDModule>
meshFittingAlgorithmMHD2DProjectionProvider("MeshFittingAlgorithmMHD2DProjection");
      
MethodCommandProvider<MeshFittingAlgorithm<MHD3DProjectionVarSet>, 
		      DataProcessingData, 
		      FiniteVolumeMHDModule>
meshFittingAlgorithmMHD3DProjectionProvider("MeshFittingAlgorithmMHD3DProjection");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
