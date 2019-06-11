#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolumeMultiFluidMHD/CoronalOutlet.hh"
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

MethodCommandProvider<CoronalOutlet, CellCenterFVMData, FiniteVolumeModule>
CoronalOutletFVMCCProvider("CoronalOutletFVMCC");

//////////////////////////////////////////////////////////////////////////////

CoronalOutlet::CoronalOutlet(const std::string& name) :
  FVMCC_BC(name)
{
}

//////////////////////////////////////////////////////////////////////////////

CoronalOutlet::~CoronalOutlet()
{
}

//////////////////////////////////////////////////////////////////////////////

void CoronalOutlet::setGhostState(GeometricEntity *const face)
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

  CFreal BxI = (*innerState)[0];
  CFreal ByI = (*innerState)[1];
  CFreal BzI = (*innerState)[2];
  CFreal VxI = (*innerState)[9];
  CFreal VyI = (*innerState)[10];
  CFreal VzI = (*innerState)[11];
  CFreal rhoI = (*innerState)[8];
  CFreal TI = (*innerState)[12];



  CFreal BrI = std::sin(thetaI)*std::cos(phiI)*BxI + std::sin(thetaI)*std::sin(phiI)*ByI + std::cos(thetaI)*BzI;
  CFreal BthetaI = std::cos(thetaI)*std::cos(phiI)*BxI + std::cos(thetaI)*std::sin(phiI)*ByI - std::sin(thetaI)*BzI;
  CFreal BphiI = -std::sin(phiI)*BxI + std::cos(phiI)*ByI;

  CFreal VrI = std::sin(thetaI)*std::cos(phiI)*VxI + std::sin(thetaI)*std::sin(phiI)*VyI + std::cos(thetaI)*VzI;
  CFreal VthetaI = std::cos(thetaI)*std::cos(phiI)*VxI + std::cos(thetaI)*std::sin(phiI)*VyI - std::sin(thetaI)*VzI;
  CFreal VphiI = -std::sin(phiI)*VxI + std::cos(phiI)*VyI;


  // Continuous BCs for r^2*Br, Btheta, r*Bphi:
  CFreal BrG = rI*rI/(rG*rG)*BrI;
  CFreal BthetaG = BthetaI;
  CFreal BphiG = (rI/rG)*BphiI;

  // Continuous BC for r^2*rho:
  CFreal rhoG = rI*rI/(rG*rG)*rhoI;

  // Continuous BCs for r^2*rho*Vr, rho*Vtheta, r*Vphi:
  CFreal VrG = rI*rI*VrI*rhoI/(rG*rG*rhoG);
  CFreal VthetaG = rhoI*VthetaI/rhoG;
  CFreal VphiG = (rI/rG)*VphiI;


  // Convert all spherical vector components back to cartesian ones:
  CFreal BxG = std::sin(thetaG)*std::cos(phiG)*BrG + std::cos(thetaG)*std::cos(phiG)*BthetaG - std::sin(phiG)*BphiG;
  CFreal ByG = std::sin(thetaG)*std::sin(phiG)*BrG + std::cos(thetaG)*std::sin(phiG)*BthetaG + std::cos(phiG)*BphiG;
  CFreal BzG = std::cos(thetaG)*BrG - std::sin(thetaG)*BthetaG;
  CFreal VxG = std::sin(thetaG)*std::cos(phiG)*VrG + std::cos(thetaG)*std::cos(phiG)*VthetaG - std::sin(phiG)*VphiG;
  CFreal VyG = std::sin(thetaG)*std::sin(phiG)*VrG + std::cos(thetaG)*std::sin(phiG)*VthetaG + std::cos(phiG)*VphiG;
  CFreal VzG = std::cos(thetaG)*VrG - std::sin(thetaG)*VthetaG;


  CFreal Bxbound = (BxI + BxG)/2.0;
  CFreal Bybound = (ByI + ByG)/2.0;
  CFreal Bzbound = (BzI + BzG)/2.0;
  CFreal Vxbound = (VxI + VxG)/2.0;
  CFreal Vybound = (VyI + VyG)/2.0;
  CFreal Vzbound = (VzI + VzG)/2.0;


  (*ghostState)[0] = BxG;
  (*ghostState)[1] = ByG;
  (*ghostState)[2] = BzG;

  // Electric field = B x V
  (*ghostState)[3] = 2.*(Bybound*Vzbound - Bzbound*Vybound) - (*innerState)[3];
  (*ghostState)[4] = 2.*(Bzbound*Vxbound - Bxbound*Vzbound) - (*innerState)[4];
  (*ghostState)[5] = 2.*(Bxbound*Vybound - Bybound*Vxbound) - (*innerState)[5];

  // According to Alejandro's setups:
  (*ghostState)[6] = -(*innerState)[6];
  (*ghostState)[7] = (*innerState)[7]; // - in Alejandro"s setup


  (*ghostState)[8] = rhoG;

  (*ghostState)[9] = VxG;
  (*ghostState)[10] = VyG;
  (*ghostState)[11] = VzG;

  // Continuous temperature
  (*ghostState)[12] = TI;


  CFLog(DEBUG_MIN, "CoronalOutlet::setGhostState() => " << (*ghostState) <<"\n");

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
