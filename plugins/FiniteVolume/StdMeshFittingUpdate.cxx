#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/StdMeshFittingUpdate.hh"
#include "FiniteVolume/FVMCC_PolyRec.hh"

#include "Framework/ComputeDummyStates.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/VolumeCalculator.hh"
#include "Framework/Node.hh"
#include "Framework/SetElementStateCoord.hh"
#include "Framework/ComputeFaceNormalsFVMCC.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdMeshFittingUpdate, CellCenterFVMData, FiniteVolumeModule> 
StdMeshFittingUpdateProvider("StdMeshFittingUpdate");

//////////////////////////////////////////////////////////////////////////////

StdMeshFittingUpdate::StdMeshFittingUpdate(const std::string& name) :
  CellCenterFVMCom(name)
{
}

//////////////////////////////////////////////////////////////////////////////

StdMeshFittingUpdate::~StdMeshFittingUpdate()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshFittingUpdate::execute()
{
  CFLog(VERBOSE, "StdMeshFittingUpdate::execute()\n");

  getMethodData().getGeoDataComputer()->modifyOffMeshNodes();
  // recompute normals, ghost state positions, volumes 
  getMethodData().getGeoDataComputer()->compute();
  // update the geometric weights for the high-order reconstruction
  getMethodData().getPolyReconstructor()->updateWeights();
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
