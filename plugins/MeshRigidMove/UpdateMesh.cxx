#include "MeshRigidMove/MeshRigidMove.hh"


#include "UpdateMesh.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshRigidMove {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateMesh, RigidMoveData, MeshRigidMoveModule> updateMeshProvider("UpdateMesh");

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UpdateMesh::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes  = socket_nodes.getDataHandle();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbNodes = nodes.size();

  RealVector rotationCenter(nbDim);
  RealVector rotationCenterTranslation(nbDim);
  RealVector translation(nbDim);
  RealVector coord(nbDim);
  RealVector expansionRatio(nbDim);

  CFreal dt = SubSystemStatusStack::getActive()->getDT();
  CFreal rotation = getMethodData().getRotationSpeed();
  rotationCenter = getMethodData().getRotationCenter();
  rotationCenterTranslation = getMethodData().getRotationCenterSpeed();
  translation = getMethodData().getTranslationSpeed();
  expansionRatio = getMethodData().getExpansionRatio();

  // Limit to two dimensional cases for rotations
  if (rotation != 0.) cf_assert (nbDim < 3);

  // Rotation angle during the timestep dt
  rotation *= dt;
  CFreal cosT = cos(rotation);
  CFreal sinT = sin(rotation);

  // Translation during the timestep dt
  translation *= dt;
  rotationCenterTranslation *= dt;

  // Translation of the center of rotation
  for (CFuint j = 0; j < nbDim; ++j) {
    rotationCenter[j] += rotationCenterTranslation[j];
  }

  for (CFuint i = 0; i < nbNodes; ++i) {
    for (CFuint j = 0; j < nbDim; ++j) {
      // Get the old coordinates
      coord[j] = (*nodes[i])[j];

      // Modify the coordinates: Translation
      coord[j] += translation[j];

      // Modify the coordinates: Expansion-Shrinking
      coord[j] = rotationCenter[j] + expansionRatio[j]*(coord[j]-rotationCenter[j]);
      }

    // Modify the coordinates: Rotation
    coord[0] = rotationCenter[0] + cosT*(coord[0]-rotationCenter[0]) + sinT*(coord[1]-rotationCenter[1]);
    coord[1] = rotationCenter[1] - sinT*(coord[0]-rotationCenter[0]) + cosT*(coord[1]-rotationCenter[1]);

    for (CFuint j = 0; j < nbDim; ++j) {
      // Set the new coordinates
      (*nodes[i])[j] = coord[j];
      }
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD
