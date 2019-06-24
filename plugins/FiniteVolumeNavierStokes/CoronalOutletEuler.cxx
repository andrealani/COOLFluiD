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
  //const CFreal zI = (dim == DIM_3D) ? innerState->getCoordinates()[ZZ] : 0.;  // we won't need the 2-D case at least for now
  const CFreal zI = innerState->getCoordinates()[ZZ];
  const CFreal rI = innerState->getCoordinates().norm2();
  const CFreal thetaI = std::atan2(std::sqrt(xI*xI + yI*yI),zI);
  const CFreal phiI = std::atan2(yI,xI);

  const CFreal xG = ghostState->getCoordinates()[XX];
  const CFreal yG = ghostState->getCoordinates()[YY];
  //const CFreal zG = (dim == DIM_3D) ? ghostState->getCoordinates()[ZZ] : 0.;
  const CFreal zG = ghostState->getCoordinates()[ZZ];
  const CFreal rG = ghostState->getCoordinates().norm2();
  const CFreal thetaG = std::atan2(std::sqrt(xG*xG + yG*yG),zG);
  const CFreal phiG = std::atan2(yG,xG);

  CFreal VxI = (*innerState)[1];
  CFreal VyI = (*innerState)[2];
  CFreal VzI = (*innerState)[3];
  CFreal rhoI = (*innerState)[0];
  CFreal TI = (*innerState)[4];
  CFreal VrI = std::sin(thetaI)*std::cos(phiI)*VxI + std::sin(thetaI)*std::sin(phiI)*VyI + std::cos(thetaI)*VzI;
  CFreal VthetaI = std::cos(thetaI)*std::cos(phiI)*VxI + std::cos(thetaI)*std::sin(phiI)*VyI - std::sin(thetaI)*VzI;
  CFreal VphiI = -std::sin(phiI)*VxI + std::cos(phiI)*VyI;
  
  // Continuous BC for r^2*rho:
  CFreal rhoG = rI*rI/(rG*rG)*rhoI;

  // Continuous BCs for r^2*rho*Vr, rho*Vtheta, r*Vphi:
  CFreal VrG = rI*rI*VrI*rhoI/(rG*rG*rhoG);
  CFreal VthetaG = rhoI*VthetaI/rhoG;
  CFreal VphiG = (rI/rG)*VphiI;
  
  CFreal VxG = std::sin(thetaG)*std::cos(phiG)*VrG + std::cos(thetaG)*std::cos(phiG)*VthetaG - std::sin(phiG)*VphiG;
  CFreal VyG = std::sin(thetaG)*std::sin(phiG)*VrG + std::cos(thetaG)*std::sin(phiG)*VthetaG + std::cos(phiG)*VphiG;
  CFreal VzG = std::cos(thetaG)*VrG - std::sin(thetaG)*VthetaG;

  
  CFreal Vxbound = (VxI + VxG)/2.0;
  CFreal Vybound = (VyI + VyG)/2.0;
  CFreal Vzbound = (VzI + VzG)/2.0;

  (*ghostState)[0] = rhoG;

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
