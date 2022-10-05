#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "SuperOutletMHD3DProjection.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperOutletMHD3DProjection, CellCenterFVMData, FiniteVolumeMHDModule> 
superOutletMHD3DProjectionFVMCCProvider("SuperOutletMHD3DProjectionFVMCC");
    
//////////////////////////////////////////////////////////////////////////////
   
void SuperOutletMHD3DProjection::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the supersonic outlet.");
   options.addConfigOption< bool >("coronalBC","Flag telling whether the coronal BC treatment is applied.");
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletMHD3DProjection::SuperOutletMHD3DProjection(const std::string& name) : 
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _refPhi = 0.;
   setParameter("refPhi",&_refPhi);

   _coronalBC = true;
   setParameter("coronalBC",&_coronalBC);
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletMHD3DProjection::~SuperOutletMHD3DProjection() 
{
}

//////////////////////////////////////////////////////////////////////

void SuperOutletMHD3DProjection::setup()
{
 FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletMHD3DProjection::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletMHD3DProjection::setGhostState(GeometricEntity *const face)
 {
   State *const innerState = face->getState(0);
   State *const ghostState = face->getState(1);
   // watch out NOT to use the operator=, because in that case 
   // the overloaded version operator=(State) would be used =>
   // also the coordinates (Node) would be set equal!!! 
   
   // set the physical data starting from the inner state
   _varSet->computePhysicalData(*innerState, _dataInnerState);

   _dataGhostState = _dataInnerState;

   if (_coronalBC) {

  CFreal RSun = 6.9551e8; // m
  //CFreal RSS = 1.4953465e10; // m
  CFreal density_code = 1.67e-13;
  CFreal velocity_code = 2.2e-4/std::sqrt(1.256e-6*1.67e-13);
  CFreal pressure_code = std::pow(2.2e-4,2)/1.256e-6;
  CFreal B_code = 2.2e-4;

  const CFreal xI = innerState->getCoordinates()[XX]*RSun;
  const CFreal xI_dimless = innerState->getCoordinates()[XX];
  //std::cout << "xI in phys. units = " << xI << endl;
  const CFreal yI = innerState->getCoordinates()[YY]*RSun;
  const CFreal yI_dimless = innerState->getCoordinates()[YY];
  const CFreal zI = innerState->getCoordinates()[ZZ]*RSun;
  const CFreal zI_dimless = innerState->getCoordinates()[ZZ];
  const CFreal rI = innerState->getCoordinates().norm2()*RSun;
  const CFreal rI_dimless = innerState->getCoordinates().norm2();
  //const CFreal thetaI = std::atan2(std::sqrt(xI*xI + yI*yI),zI);
  //const CFreal phiI = std::atan2(yI,xI);
  const CFreal rhoI = std::sqrt(xI*xI + yI*yI);
  const CFreal rhoI_dimless = std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless);

  const CFreal xG = ghostState->getCoordinates()[XX]*RSun;
  const CFreal xG_dimless = ghostState->getCoordinates()[XX];
  const CFreal yG = ghostState->getCoordinates()[YY]*RSun;
  const CFreal yG_dimless = ghostState->getCoordinates()[YY];
  const CFreal zG = ghostState->getCoordinates()[ZZ]*RSun;
  const CFreal zG_dimless = ghostState->getCoordinates()[ZZ];
  const CFreal rG = ghostState->getCoordinates().norm2()*RSun;
  const CFreal rG_dimless = ghostState->getCoordinates().norm2();
  //const CFreal thetaG = std::atan2(std::sqrt(xG*xG + yG*yG),zG);
  //const CFreal phiG = std::atan2(yG,xG);
  const CFreal rhoG = std::sqrt(xG*xG + yG*yG);
  const CFreal rhoG_dimless = std::sqrt(xG_dimless*xG_dimless + yG_dimless*yG_dimless);

  CFreal densityI = (*innerState)[0]*density_code;
  CFreal densityI_dimless = (*innerState)[0];
  //std::cout << "density in phys. units = " << densityI << endl;
  CFreal VxI = (*innerState)[1]*velocity_code;
  CFreal VxI_dimless = (*innerState)[1];
  //std::cout << "Vx in phys. units = " << VxI << endl;
  CFreal VyI = (*innerState)[2]*velocity_code;
  CFreal VyI_dimless = (*innerState)[2];
  CFreal VzI = (*innerState)[3]*velocity_code;
  CFreal VzI_dimless = (*innerState)[3];
  CFreal BxI = (*innerState)[4]*B_code;
  CFreal BxI_dimless = (*innerState)[4];
  //std::cout << "Bx in phys. units = " << BxI << endl;
  CFreal ByI = (*innerState)[5]*B_code;
  CFreal ByI_dimless = (*innerState)[5];
  CFreal BzI = (*innerState)[6]*B_code;
  CFreal BzI_dimless = (*innerState)[6];
  CFreal PrI = (*innerState)[7]*pressure_code;
  CFreal PrI_dimless = (*innerState)[7];
  //std::cout << "P in phys. units = " << PrI << endl;



  CFreal BrI_dimless = xI_dimless/rI_dimless*BxI_dimless + yI_dimless/rI_dimless*ByI_dimless + zI_dimless/rI_dimless*BzI_dimless;
  CFreal BthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*BxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*ByI_dimless - rhoI_dimless/rI_dimless*BzI_dimless;
  CFreal BphiI_dimless = -yI_dimless/rhoI_dimless*BxI_dimless + xI_dimless/rhoI_dimless*ByI_dimless;

  CFreal VrI_dimless = xI_dimless/rI_dimless*VxI_dimless + yI_dimless/rI_dimless*VyI_dimless + zI_dimless/rI_dimless*VzI_dimless;
  CFreal VthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VyI_dimless - rhoI_dimless/rI_dimless*VzI_dimless;
  CFreal VphiI_dimless = -yI_dimless/rhoI_dimless*VxI_dimless + xI_dimless/rhoI_dimless*VyI_dimless;




  // Continuous BCs for r^2*Br, Btheta, Bphi:
  CFreal BrG_dimless = rI_dimless*rI_dimless/(rG_dimless*rG_dimless)*BrI_dimless;
  CFreal BthetaG_dimless = BthetaI_dimless;
  //Jens:
  CFreal BphiG_dimless = (rI_dimless/rG_dimless)*BphiI_dimless;
  //CFreal BphiG_dimless = BphiI_dimless;
  
  CFreal densityG_dimless = rI_dimless*rI_dimless*rI_dimless*densityI_dimless/(rG_dimless*rG_dimless*rG_dimless);
  (*ghostState)[0] = densityG_dimless;

  // Continuous BCs for r^2*density*Vr, density*Vtheta, r*Vphi:
  CFreal VrG_dimless = rI_dimless*rI_dimless*VrI_dimless*densityI_dimless/(rG_dimless*rG_dimless*densityG_dimless);
  //CFreal VrG_dimless = rI_dimless*rI_dimless*VrI_dimless/(rG_dimless*rG_dimless);
  CFreal VthetaG_dimless = densityI_dimless*VthetaI_dimless/densityG_dimless;
  CFreal VphiG_dimless = (rI_dimless/rG_dimless)*VphiI_dimless;

  // UPDATE: BARBARA'S VELOCITY OUTLET CONDITIONS:
  //CFreal VrG_dimless = VrI_dimless;
  //CFreal VrG_dimless = rI_dimless*rI_dimless*VrI_dimless/(rG_dimless*rG_dimless);
  //CFreal VthetaG_dimless = VthetaI_dimless;
  //CFreal VphiG_dimless = VphiI_dimless;


  // Convert all spherical vector components back to cartesian ones:
  CFreal BxG_dimless = xG_dimless/rG_dimless*BrG_dimless + xG_dimless*zG_dimless/(rhoG_dimless*rG_dimless)*BthetaG_dimless - yG_dimless/rhoG_dimless*BphiG_dimless;
  CFreal ByG_dimless = yG_dimless/rG_dimless*BrG_dimless + yG_dimless*zG_dimless/(rhoG_dimless*rG_dimless)*BthetaG_dimless + xG_dimless/rhoG_dimless*BphiG_dimless;
  CFreal BzG_dimless = zG_dimless/rG_dimless*BrG_dimless - rhoG_dimless/rG_dimless*BthetaG_dimless;
  CFreal VxG_dimless = xG_dimless/rG_dimless*VrG_dimless + xG_dimless*zG_dimless/(rhoG_dimless*rG_dimless)*VthetaG_dimless - yG_dimless/rhoG_dimless*VphiG_dimless;
  CFreal VyG_dimless = yG_dimless/rG_dimless*VrG_dimless + yG_dimless*zG_dimless/(rhoG_dimless*rG_dimless)*VthetaG_dimless + xG_dimless/rhoG_dimless*VphiG_dimless;
  CFreal VzG_dimless = zG_dimless/rG_dimless*VrG_dimless - rhoG_dimless/rG_dimless*VthetaG_dimless;


/*
  CFreal Bxbound = (BxI + BxG)/2.0;
  CFreal Bybound = (ByI + ByG)/2.0;
  CFreal Bzbound = (BzI + BzG)/2.0;
  CFreal Vxbound = (VxI + VxG)/2.0;
  CFreal Vybound = (VyI + VyG)/2.0;
  CFreal Vzbound = (VzI + VzG)/2.0;

  // Boundary values for the E-field:
  CFreal ExB = Bybound*Vzbound - Bzbound*Vybound;
  CFreal EyB = Bzbound*Vxbound - Bxbound*Vzbound;
  CFreal EzB = Bxbound*Vybound - Bybound*Vxbound;
*/

  (*ghostState)[4] = BxG_dimless;
  (*ghostState)[5] = ByG_dimless;
  (*ghostState)[6] = BzG_dimless;



  // density*r = continuous => densityG = densityI*rI/rG
  // Continuous BC for r^2*rho: TYPO IN SKRALAN's PAPER
  //CFreal densityG = rI*rI/(rG*rG)*densityI;


  (*ghostState)[1] = VxG_dimless;
  (*ghostState)[2] = VyG_dimless;
  (*ghostState)[3] = VzG_dimless;

  
  (*ghostState)[7] = (*innerState)[7]*rI_dimless*rI_dimless*rI_dimless/(rG_dimless*rG_dimless*rG_dimless);

   //
   //}

  (*ghostState)[8] = -(*innerState)[8];
  }
  else {
   // there are two possible boundary conditions for phi

   // 1) a reference value should be imposed for phi
   _dataGhostState[MHDProjectionTerm::PHI] = _refPhi;

   // 2)
   //_dataGhostState[MHDProjectionTerm::PHI] = -_dataInnerState[MHDProjectionTerm::PHI];

   _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
  }
 }

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
