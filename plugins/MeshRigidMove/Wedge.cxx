#include "MeshRigidMove/MeshRigidMove.hh"


#include "Wedge.hh"
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

MethodCommandProvider<Wedge, RigidMoveData, MeshRigidMoveModule> wedgeProvider("Wedge");

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
Wedge::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Wedge::execute()
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

  CFreal dt = SubSystemStatusStack::getActive()->getDTDim();

  for (CFuint i = 0; i < nbNodes; ++i) {
    for (CFuint j = 0; j < nbDim; ++j) {
      // Get the old coordinates
      coord[j] = (*nodes[i])[j];
      }

    /// Test modification (Wedge testcase)
    CFreal limit = 0.5;
    CFreal ratio = 5.;
    if(SubSystemStatusStack::getActive()->getCurrentTimeDim() > 0.7){
    if (coord[0] > limit-(coord[1]/3+(100*ratio*dt))){
      coord[0] = coord[0];
      coord[1] += ((1.-coord[1])*(1/(3*(1.+coord[1])))*ratio*dt*coord[1]);
      }
    }

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
