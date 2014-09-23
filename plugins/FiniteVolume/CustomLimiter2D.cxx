#include "CustomLimiter2D.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"
#include "Common/CFLog.hh"
#include "MathTools/MathFunctions.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<CustomLimiter2D,CellCenterFVMData,
		       Limiter<CellCenterFVMData>,FiniteVolumeModule>
customLimiter2DProvider("Custom2D");

//////////////////////////////////////////////////////////////////////////////

CustomLimiter2D::CustomLimiter2D(const std::string& name) :
  CustomLimiter1D(name), 
  socket_uY("uY")
{
}

//////////////////////////////////////////////////////////////////////////////

CustomLimiter2D::~CustomLimiter2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void CustomLimiter2D::computeDeltaMin(const Node& coord, const State& state, CFuint iVar)
{
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  const CFuint stateID = state.getLocalID();
  const RealVector& stateCoord = state.getCoordinates();
  _deltaMin = uX(stateID,iVar,state.size())*(coord[XX] - stateCoord[XX]) + 
    uY(stateID,iVar,state.size())*(coord[YY] - stateCoord[YY]);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

