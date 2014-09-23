#include "FiniteVolume/FiniteVolume.hh"
#include "BDF2ALEUpdate.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<BDF2ALEUpdate, CellCenterFVMData, FiniteVolumeModule> BDF2ALEUpdateProvider("BDF2ALEUpdate");

//////////////////////////////////////////////////////////////////////////////

void BDF2ALEUpdate::updateNormalsData()
{
  CFAUTOTRACE;

  // Update the normals
  StdALEUpdate::updateNormalsData();

  // Store the normals
  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  DataHandle<CFreal> avNormals = socket_avNormals.getDataHandle();

  CFuint nbNormals = normals.size();
  CFreal temp;

  const CFreal alpha = SubSystemStatusStack::getActive()->getPreviousDT() / SubSystemStatusStack::getActive()->getDT();
  const CFreal xi = 1./(1. + alpha);

  const bool isFirstTimeStep = (SubSystemStatusStack::getActive()->getNbIter() == 1) ?
    true : false;

  if(!isFirstTimeStep){
    for (CFuint iNormal = 0; iNormal < nbNormals; ++iNormal) {
      temp = (1.+ xi) * normals[iNormal] - xi*avNormals[iNormal];
      avNormals[iNormal] = normals[iNormal];
      normals[iNormal] = temp;
    }
  }
  else{
    for (CFuint iNormal = 0; iNormal < nbNormals; ++iNormal) {
      avNormals[iNormal] = normals[iNormal];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
BDF2ALEUpdate::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result =
    StdALEUpdate::needsSockets();

  result.push_back(&socket_avNormals);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
