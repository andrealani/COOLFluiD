#include "MeshRigidMove/MeshRigidMove.hh"


#include "Piston.hh"
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

MethodCommandProvider<Piston, RigidMoveData, MeshRigidMoveModule> pistonProvider("Piston");

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
Piston::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Piston::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes  = socket_nodes.getDataHandle();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbNodes = nodes.size();
  RealVector coord(nbDim);

  const CFreal dt = SubSystemStatusStack::getActive()->getDTDim();

  // Define accelaration derivative
  const CFreal accelerationDeriv = 0.2;
  // Maximum grid X coordinate
  const CFreal maximumX = 5.;

  // Compute acceleration, speed, displacement of piston
  const CFreal acceleration = _oldAcceleration + accelerationDeriv*dt;
  const CFreal speed = _oldSpeed + acceleration*dt;
  const CFreal displacement = speed*dt;
CFout << "Piston speed: " <<speed <<"\n";
  for (CFuint i = 0; i < nbNodes; ++i) {
    for (CFuint j = 0; j < nbDim; ++j) {
      // Get the old coordinates
      coord[j] = (*nodes[i])[j];
    }

    /// Test modification (piston testcase)
    coord[0] += displacement * (maximumX-coord[0])/(maximumX-_totalDisplacement);
    for (CFuint j = 0; j < nbDim; ++j) {
      // Set the new coordinates
      (*nodes[i])[j] = coord[j];
    }

  }

  // Compute the new position of piston
  _totalDisplacement += displacement;
  _oldAcceleration = acceleration;
  _oldSpeed = speed;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD
