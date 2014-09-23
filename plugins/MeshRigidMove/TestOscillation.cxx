#include "MeshRigidMove/MeshRigidMove.hh"


#include "TestOscillation.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/Node.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshRigidMove {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TestOscillation, RigidMoveData, MeshRigidMoveModule> TestOscillationProvider("TestOscillation");

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TestOscillation::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

/// @todo Declare the private data !!!!!

void TestOscillation::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes  = socket_nodes.getDataHandle();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbNodes = nodes.size();

  RealVector rotationCenter(nbDim);
  RealVector coord(nbDim);

  rotationCenter = getMethodData().getRotationCenter();

  // Limit to two dimensional cases for rotations
  cf_assert (nbDim < 3);

  // Compute Rotation Angle at Current Time
  CFreal time = SubSystemStatusStack::getActive()->getCurrentTimeDim();
  CFreal rotation = (0.016 + 2.51*(sin(2.*0.0814*time)))*100.;

  // Rotation Angle = CurrentAlpha - PastAlpha
  rotation -= _currentAlpha;


  for (CFuint i = 0; i < nbNodes; ++i) {

    for (CFuint j = 0; j < nbDim; ++j) {
      // Get the old coordinates
      coord[j] = (*nodes[i])[j];
    }

    CFreal distanceToCenter = sqrt((coord[0]-rotationCenter[0])*(coord[0]-rotationCenter[0]) + (coord[1]-rotationCenter[1])*(coord[1]-rotationCenter[1]));
distanceToCenter = min(distanceToCenter/.5, 1.);

    CFreal backrotation = rotation;
    rotation *= (1.- distanceToCenter);

    // Rotation angle
    CFreal cosT = cos(rotation*MathTools::MathConsts::CFrealPi()/180.);
    CFreal sinT = sin(rotation*MathTools::MathConsts::CFrealPi()/180.);

    // Modify the coordinates: Rotation
    coord[0] = rotationCenter[0] + cosT*(coord[0]-rotationCenter[0]) + sinT*(coord[1]-rotationCenter[1]);
    coord[1] = rotationCenter[1] - sinT*(coord[0]-rotationCenter[0]) + cosT*(coord[1]-rotationCenter[1]);

    rotation = backrotation;
    for (CFuint j = 0; j < nbDim; ++j) {
      // Set the new coordinates
      (*nodes[i])[j] = coord[j];
    }
  }

  // Compute the new position of TestOscillation
  _currentAlpha += rotation;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD
