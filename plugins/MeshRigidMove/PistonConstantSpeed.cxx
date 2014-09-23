#include "MeshRigidMove/MeshRigidMove.hh"


#include "PistonConstantSpeed.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
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

MethodCommandProvider<PistonConstantSpeed, RigidMoveData, MeshRigidMoveModule> pistonConstantSpeedProvider("PistonConstantSpeed");

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
PistonConstantSpeed::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

/// @todo Declare the private data !!!!!

void PistonConstantSpeed::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes  = socket_nodes.getDataHandle();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbNodes = nodes.size();
  RealVector coord(nbDim);

  CFreal dt = SubSystemStatusStack::getActive()->getDTDim();

  // Maximum grid X coordinate
  CFreal maximumX = 5.;

  // Compute acceleration, speed, displacement of PistonConstantSpeed
  CFreal speed = 0.8;
  CFreal displacement = speed*dt;

  for (CFuint i = 0; i < nbNodes; ++i) {
    for (CFuint j = 0; j < nbDim; ++j) {
      // Get the old coordinates
      coord[j] = (*nodes[i])[j];
      }

    /// PistonConstantSpeed testcase
    coord[0] += displacement * (maximumX-coord[0])/(maximumX-_totalDisplacement);

    for (CFuint j = 0; j < nbDim; ++j) {
      // Set the new coordinates
      (*nodes[i])[j] = coord[j];
      }
  }

  // Compute the new position of PistonConstantSpeed
  _totalDisplacement += displacement;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD
