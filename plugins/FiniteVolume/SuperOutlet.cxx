#include "FiniteVolume/FiniteVolume.hh"
#include "SuperOutlet.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperOutlet, CellCenterFVMData, FiniteVolumeModule> superOutletFVMCCProvider("SuperOutletFVMCC");

//////////////////////////////////////////////////////////////////////////////

SuperOutlet::SuperOutlet(const std::string& name) :
  FVMCC_BC(name)
{
}

//////////////////////////////////////////////////////////////////////////////

SuperOutlet::~SuperOutlet()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutlet::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  // watch out NOT to use the operator=, because in that case
  // the overloaded version operator=(State) would be used =>
  // also the coordinates (Node) would be set equal!!!
  ghostState->copyData(*innerState);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
