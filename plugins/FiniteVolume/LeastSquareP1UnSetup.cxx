#include "FiniteVolume/FiniteVolume.hh"
#include "LeastSquareP1UnSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LeastSquareP1UnSetup, CellCenterFVMData, FiniteVolumeModule> leastSquareP1UnSetupProvider("LeastSquareP1UnSetup");

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1UnSetup::execute()
{
  CFAUTOTRACE;

  StdUnSetup::execute();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
