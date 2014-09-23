#include "CustomLimiter3D.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"
#include "Common/CFLog.hh"
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

MethodStrategyProvider<CustomLimiter3D,CellCenterFVMData,
		       Limiter<CellCenterFVMData>,FiniteVolumeModule>
customLimiter3DProvider("Custom3D");

//////////////////////////////////////////////////////////////////////////////

CustomLimiter3D::CustomLimiter3D(const std::string& name) :
  CustomLimiter2D(name),
  socket_uZ("uZ")
{
}

//////////////////////////////////////////////////////////////////////////////

CustomLimiter3D::~CustomLimiter3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void CustomLimiter3D::computeDeltaMin(const Node& coord, 
				      const State& state, 
				      CFuint iVar)
{
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> uZ = socket_uZ.getDataHandle();
  
  const CFuint stateID = state.getLocalID();
  const RealVector& stateCoord = state.getCoordinates();
  _deltaMin = uX(stateID,iVar,state.size())*(coord[XX] - stateCoord[XX]) + 
    uY(stateID,iVar,state.size())*(coord[YY] - stateCoord[YY]) + 
    uZ(stateID,iVar,state.size())*(coord[ZZ] - stateCoord[ZZ]);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

