#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/CoronalOutletEuler.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CoronalOutletEuler, CellCenterFVMData, FiniteVolumeNavierStokesModule>
CoronalOutletEulerFVMCCProvider("CoronalOutletEulerFVMCC");

//////////////////////////////////////////////////////////////////////////////

CoronalOutletEuler::CoronalOutletEuler(const std::string& name) :
  FVMCC_BC(name)
{
}

//////////////////////////////////////////////////////////////////////////////

CoronalOutletEuler::~CoronalOutletEuler()
{
}

//////////////////////////////////////////////////////////////////////////////

void CoronalOutletEuler::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  // watch out NOT to use the operator=, because in that case
  // the overloaded version operator=(State) would be used =>
  // also the coordinates (Node) would be set equal!!!
  ghostState->copyData(*innerState);

  //const CFuint dim = _bCoord.size();
  const CFuint dim = DIM_3D;   // only preliminary as _bCoord not available!
  cf_assert(dim == DIM_3D);

  const CFreal xI = innerState->getCoordinates()[XX];
  const CFreal yI = innerState->getCoordinates()[YY];
  const CFreal zI = innerState->getCoordinates()[ZZ];
  const CFreal rI = innerState->getCoordinates().norm2();
  const CFreal rhoI = std::sqrt(xI*xI + yI*yI);
  const CFreal thetaI = std::atan2(std::sqrt(xI*xI + yI*yI),zI);
  const CFreal phiI = std::atan2(yI,xI);

  const CFreal xG = ghostState->getCoordinates()[XX];
  const CFreal yG = ghostState->getCoordinates()[YY];
  const CFreal zG = ghostState->getCoordinates()[ZZ];
  const CFreal rG = ghostState->getCoordinates().norm2();
  const CFreal thetaG = std::atan2(std::sqrt(xG*xG + yG*yG),zG);
  const CFreal phiG = std::atan2(yG,xG);
  const CFreal rhoG = std::sqrt(xG*xG + yG*yG);

  CFreal VxI = (*innerState)[1];
  CFreal VyI = (*innerState)[2];
  CFreal VzI = (*innerState)[3];
  CFreal densityI = (*innerState)[0];
  CFreal TI = (*innerState)[4];
  CFreal VrI = xI/rI*VxI + yI/rI*VyI + zI/rI*VzI;
  CFreal VthetaI = xI*zI/(rhoI*rI)*VxI + yI*zI/(rhoI*rI)*VyI - rhoI/rI*VzI;
  CFreal VphiI = -yI/rhoI*VxI + xI/rhoI*VyI;
  
  // Continuous BC for r^2*rho:
  CFreal densityG = rI*rI/(rG*rG)*densityI;

  // Continuous BCs for r^2*rho*Vr, rho*Vtheta, r*Vphi:
  CFreal VrG = rI*rI*VrI*densityI/(rG*rG*densityG);
  //CFreal VrG = rI*rI*VrI/(rG*rG);
  CFreal VthetaG = VthetaI; //rhoI*VthetaI/densityG;
  CFreal VphiG = (rI/rG)*VphiI;
  
 CFreal VxG = xG/rG*VrG - yG/rhoG*VphiG + xG*zG/(rhoG*rG)*VthetaG;
 CFreal VyG = yG/rG*VrG + xG/rhoG*VphiG + yG*zG/(rhoG*rG)*VthetaG;
 CFreal VzG = zG/rG*VrG - rhoG/rG*VthetaG;

  
  CFreal Vxbound = (VxI + VxG)/2.0;
  CFreal Vybound = (VyI + VyG)/2.0;
  CFreal Vzbound = (VzI + VzG)/2.0;

  (*ghostState)[0] = densityG;

  (*ghostState)[1] = VxG;
  (*ghostState)[2] = VyG;
  (*ghostState)[3] = VzG;
  
  // Continuous temperature
  (*ghostState)[4] = TI;
  

  CFLog(DEBUG_MIN, "CoronalOutletEuler::setGhostState() => " << (*ghostState) <<"\n");

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
