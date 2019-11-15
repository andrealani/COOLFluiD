#include "FiniteVolume/FiniteVolume.hh"
#include "SuperInletProjection.hh"
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

MethodCommandProvider<SuperInletProjection, CellCenterFVMData, FiniteVolumeModule> 
superInletProjectionFVMCCProvider("SuperInletProjectionFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjection::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the supersonic inlet.");
  options.addConfigOption< vector<CFuint> >("ProjectionIDs","IDs corresponding to projection fields.");
  options.addConfigOption< CFint >("inletCoronalBC","Switch to coronal inlet BC.");
  options.addConfigOption< CFint >("Phi_divB_zero","Set Phi to zero at the ghost state.");
  options.addConfigOption< CFint >("Phi_divB_extrapolated","Copy Phi from the closest inner state to the ghost cell.");
}
      
//////////////////////////////////////////////////////////////////////////////

SuperInletProjection::SuperInletProjection(const std::string& name) :
  SuperInlet(name)
{
  addConfigOptionsTo(this);

  _refPhi = 0.;
  setParameter("refPhi",&_refPhi);
  
  setParameter("inletCoronalBC",&_inletCoronalBC);
  setParameter("Phi_divB_zero",&_Phi_divB_zero);
  setParameter("Phi_divB_extrapolated",&_Phi_divB_extrapolated);
  
  _projectionIDs = vector<CFuint>();
  setParameter("ProjectionIDs",&_projectionIDs);
}

//////////////////////////////////////////////////////////////////////////////

SuperInletProjection::~SuperInletProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjection::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  if (_inletCoronalBC==1) {

  CFreal latG = 0.;
  CFreal RSun = 6.9551e8; // m
  CFreal RSS = 1.4953465e10; // m
  CFreal B0dip = 0.00022; // Tesla
  CFreal latI = 0.;
  const CFreal PI = MathTools::MathConsts::CFrealPi();
  const CFreal xG = ghostState->getCoordinates()[XX]*RSun;
  const CFreal xG_dimless = ghostState->getCoordinates()[XX];
  const CFreal yG = ghostState->getCoordinates()[YY]*RSun;
  const CFreal yG_dimless = ghostState->getCoordinates()[YY];
  const CFreal zG = ghostState->getCoordinates()[ZZ]*RSun;
  const CFreal zG_dimless = ghostState->getCoordinates()[ZZ];
  const CFreal xI = innerState->getCoordinates()[XX]*RSun;
  const CFreal xI_dimless = innerState->getCoordinates()[XX];
  const CFreal yI = innerState->getCoordinates()[YY]*RSun;
  const CFreal yI_dimless = innerState->getCoordinates()[YY];
  const CFreal zI = innerState->getCoordinates()[ZZ]*RSun;
  const CFreal zI_dimless = innerState->getCoordinates()[ZZ];
  const CFreal xBoundary = (xG + xI)/2.0;
  const CFreal yBoundary = (yG + yI)/2.0;
  const CFreal zBoundary = (zG + zI)/2.0;
  const CFreal rBoundary = std::sqrt(xBoundary*xBoundary + yBoundary*yBoundary + zBoundary*zBoundary);
  const CFreal rhoBoundary = std::sqrt(xBoundary*xBoundary + yBoundary*yBoundary);
  const CFreal thetaBoundary = std::atan2(std::sqrt(xBoundary*xBoundary + yBoundary*yBoundary),zBoundary);
  const CFreal phiBoundary = std::atan2(yBoundary,xBoundary);

/*

  // Convert boundary values provided in the CFcase file into spherical coordinates
  CFreal Vx_user = (*_dimState)[1];
  CFreal Vy_user = (*_dimState)[2];
  CFreal Vz_user = (*_dimState)[3];
  CFreal Vr_user = xBoundary/rBoundary*Vx_user + yBoundary/rBoundary*Vy_user + zBoundary/rBoundary*Vz_user;
  CFreal Vtheta_user = xBoundary*zBoundary/(rhoBoundary*rBoundary)*Vx_user + yBoundary*zBoundary/(rhoBoundary*rBoundary)*Vy_user - rhoBoundary/rBoundary*Vz_user;
  CFreal Vphi_user = -yBoundary/rhoBoundary*Vx_user + xBoundary/rhoBoundary*Vy_user;

*/


  const CFreal rG = std::sqrt(xG*xG + yG*yG + zG*zG);
  const CFreal rG_dimless = std::sqrt(xG_dimless*xG_dimless + yG_dimless*yG_dimless + zG_dimless*zG_dimless);
  const CFreal rhoG = std::sqrt(xG*xG + yG*yG);
  const CFreal rhoG_dimless = std::sqrt(xG_dimless*xG_dimless + yG_dimless*yG_dimless);
  const CFreal thetaG = std::atan2(std::sqrt(xG*xG + yG*yG),zG);
  const CFreal phiG = std::atan2(yG,xG);
  const CFreal rI = std::sqrt(xI*xI + yI*yI + zI*zI);
  const CFreal rI_dimless = std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless + zI_dimless*zI_dimless);
  const CFreal thetaI = std::atan2(std::sqrt(xI*xI + yI*yI),zI);
  const CFreal rhoI = std::sqrt(xI*xI + yI*yI);
  const CFreal rhoI_dimless = std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless);
  const CFreal phiI = std::atan2(yI,xI);

  if (thetaG > -PI && thetaG < -PI*0.5) {
     latG = std::abs(thetaG) - PI*0.5;
  } else if (thetaG > -PI*0.5 && thetaG < 0) {
     latG = PI*0.5 - std::abs(thetaG);
  } else if (thetaG > 0. && thetaG < PI*0.5) {
     latG = PI*0.5 - thetaG;
  } else if (thetaG > PI*0.5 && thetaG < PI) {
     latG = thetaG - PI*0.5;
  } else {
    std::cout << "Error: value of theta for the point in question outside expected range" << endl;
  }
  
  // DENSITY kept constant (Dirichlet condition)

  // For the electric field, Phi & Psi
  // Continuous BC for the electric field -- not working
  //(*ghostState)[3] = 2.*(*_dimState)[3] - (*innerState)[3]; //Ex
  //(*ghostState)[4] = 2.*(*_dimState)[4] - (*innerState)[4]; //Ey
  //(*ghostState)[5] = 2.*(*_dimState)[4] - (*innerState)[4]; //Ez
  
  // Zero electric field at the boundary:
  //(*ghostState)[3] = - (*innerState)[3]; //Ex
  //(*ghostState)[4] = - (*innerState)[4]; //Ey
  //(*ghostState)[5] = - (*innerState)[4]; //Ez
  
  // More consistently: Zero electric field E = - V x B at the boundary:
  // Vx = 9 Vy = 10 Vz = 11 Bx = 0 By = 1 Bz = 2
  //(*ghostState)[3] =2.*((*_dimState)[11]*(*_dimState)[1] - (*_dimState)[2]*(*_dimState)[10]) - (*innerState)[3];
  //(*ghostState)[4] = 2.*((*_dimState)[9]*(*_dimState)[2] - (*_dimState)[11]*(*_dimState)[0]) - (*innerState)[4];
  //(*ghostState)[5] = 2.*((*_dimState)[10]*(*_dimState)[0] - (*_dimState)[9]*(*_dimState)[1]) - (*innerState)[5];

  //(*ghostState)[3] = - (*innerState)[3]; //2.*((*_dimState)[11]*(*_dimState)[1] - (*_dimState)[2]*(*_dimState)[10]) - (*innerState)[3];
  //(*ghostState)[4] = - (*innerState)[4];  //2.*((*_dimState)[9]*(*_dimState)[2] - (*_dimState)[11]*(*_dimState)[0]) - (*innerState)[4];
  //(*ghostState)[5] = - (*innerState)[5];   //2.*((*_dimState)[10]*(*_dimState)[0] - (*_dimState)[9]*(*_dimState)[1]) - (*innerState)[5];
  
  // - VxB BC but V and B are from case file so give zero anyway
  //(*ghostState)[3] = 2.*((*_dimState)[11]*(*_dimState)[1] - (*_dimState)[2]*(*_dimState)[10]) - (*innerState)[3]; //Ex
  //(*ghostState)[4] = 2.*((*_dimState)[9]*(*_dimState)[2] - (*_dimState)[11]*(*_dimState)[0]) - (*innerState)[4]; //Ey
  //(*ghostState)[5] = 2.*((*_dimState)[10]*(*_dimState)[0] - (*_dimState)[9]*(*_dimState)[1]) - (*innerState)[5]; //Ex

  // Perfectly conducting wall BC
  //(*ghostState)[3] = - (*innerState)[3] + 2.0*((*innerState)[3]*std::sin(thetaI)*std::cos(phiI) + (*innerState)[4]*std::sin(thetaI)*std::sin(phiI) + (*innerState)[5]*std::cos(thetaI))*std::sin(thetaI)*std::cos(phiI);
  //(*ghostState)[4] = - (*innerState)[4] + 2.0*((*innerState)[3]*std::sin(thetaI)*std::cos(phiI) + (*innerState)[4]*std::sin(thetaI)*std::sin(phiI) + (*innerState)[5]*std::cos(thetaI))*std::sin(thetaI)*std::sin(phiI);
  //(*ghostState)[5] = - (*innerState)[5] + 2.0*((*innerState)[3]*std::sin(thetaI)*std::cos(phiI) + (*innerState)[4]*std::sin(thetaI)*std::sin(phiI) + (*innerState)[5]*std::cos(thetaI))*std::cos(thetaI);
  //(*ghostState)[3] = - (*innerState)[3];
  //(*ghostState)[4] = - (*innerState)[4];
  //(*ghostState)[5] = - (*innerState)[5];

  // Continuous condition for Phi, Psi
  //(*ghostState)[6] = (*innerState)[6];
  //(*ghostState)[7] = (*innerState)[7];



  // divB cleaning Phi BC:
  if (_Phi_divB_zero==1) {
    (*ghostState)[8] = -(*innerState)[8];
  }
  else if (_Phi_divB_extrapolated==1) {
    (*ghostState)[8] = (*innerState)[8];
  }
  else {
    std::cout << "ERROR: Phi not specified at inlet boundary!" << endl;
  }





//-(*innerState)[6]; // a + in Alejandro's case // Also zero in Dana's case





  //(*ghostState)[7] = -(*innerState)[7]; // - in Alejandro's case
  CFreal densityBoundary_code_units = 1.0; // in code units
  //CFreal densityG = 2.*(*_dimState)[0] - (*innerState)[0];
  CFreal densityG = 2.*densityBoundary_code_units - (*innerState)[0];
  (*ghostState)[0] = densityG;
  //std::cout << "rhoB_cu = " << densityBoundary_code_units << endl;
  //std::cout << "rhoI_cu = " << (*innerState)[0] << endl;
  //std::cout << "ghoststate rho = " << (*ghostState)[8] << endl;
  
  // VELOCITY: This includes differential rotation but the first version should not
  CFreal A = 14.713;   // values from Wikipedia; find a better and more accurate source...
  CFreal B = -2.396;
  CFreal C = -1.787;
  
  CFreal Omega_deg_day = A + B*std::sin(latG)*std::sin(latG) + C*std::sin(latG)*std::sin(latG)*std::sin(latG)*std::sin(latG);
  // Convert OmegaG from deg/day to rad/s
  // deg/day = (PI/180.*rad)/(24*60*60*sec)
  //OmegaG = OmegaG*(PI/180)/(24*60*60);
  //OmegaI = OmegaI*(PI/180)/(24*60*60);
  CFreal Omega = Omega_deg_day*PI/(180*24*60*60);   // Now in rad/sec


  //CFreal VrG = 0.;
  //CFreal VrI = 0.;
  //CFreal VthetaG = 0.;
  //CFreal VthetaI = 0.;
  //CFreal VxG = 0.;
  //CFreal VyG = 0.;
  //CFreal VzG = 0.;




//  CFreal Bxds = (*_dimState)[0];
// //cout << "Bxds = " << Bxds << endl;
//  CFreal Byds = (*_dimState)[1];
// //cout << "Byds = " << Byds << endl;
//  CFreal Bzds = (*_dimState)[2];
// //cout << "Bzds = " << Bzds << endl;

//  CFreal Brds = xBoundary/rBoundary*Bxds + yBoundary/rBoundary*Byds + zBoundary/rBoundary*Bzds;
// //cout << "Brds = " << Brds << endl;
//  CFreal Bthetads = xBoundary*zBoundary/(rhoBoundary*rBoundary)*Bxds + yBoundary*zBoundary/(rhoBoundary*rBoundary)*Byds - rhoBoundary/rBoundary*Bzds;
// //cout << "Bthetads = " << Bthetads << endl;
//  CFreal Bphids = -yBoundary/rhoBoundary*Bxds + xBoundary/rhoBoundary*Byds;
// //cout << "Bphids = " << Bphids << endl;

  CFreal BxI = (*innerState)[4]*2.2e-4;
 //cout << "BxI = " << BxI << endl;
  CFreal ByI = (*innerState)[5]*2.2e-4;
 //cout << "ByI = " << ByI << endl;
  CFreal BzI = (*innerState)[6]*2.2e-4;
 //cout << "BzI = " << BzI << endl;
  CFreal BrI = xI/rI*BxI + yI/rI*ByI + zI/rI*BzI;
 //cout << "BrI = " << BrI << endl;
  CFreal BthetaI = xI*zI/(rhoI*rI)*BxI + yI*zI/(rhoI*rI)*ByI - rhoI/rI*BzI;
 //cout << "BthetaI = " << BthetaI << endl;
  CFreal BphiI = -yI/rhoI*BxI + xI/rhoI*ByI;
 //cout << "BphiI = " << BphiI << endl;



  // VphiG = 2*VphiB - VphiI = -VphiI; // only preliminary before taking differential rotation into account; this value is independent of the zone
  //CFreal densityI = (*innerState)[0]*1.67e-13;
  CFreal VxI = (*innerState)[1]*2.2e-4/sqrt(1.256e-6*1.67e-13);
  CFreal VxI_dimless = (*innerState)[1];
  CFreal VyI = (*innerState)[2]*2.2e-4/sqrt(1.256e-6*1.67e-13);
  CFreal VyI_dimless = (*innerState)[2];
  CFreal VzI = (*innerState)[3]*2.2e-4/sqrt(1.256e-6*1.67e-13);
  CFreal VzI_dimless = (*innerState)[3];
  CFreal VphiI = -yI/rhoI*VxI + xI/rhoI*VyI;
  //CFreal VphiBoundary = Omega*rBoundary*std::sin(thetaBoundary);
  //CFreal VphiG = VphiI;   // extrapolated according to Jens

//std::cout << "Vphi = " << VphiG << endl;

//  if (latG  > 0.3927) {  // this value corresponds to 45/2 deg in rad
    //cout << "in the pole region" << endl;
    CFreal VrI = xI/rI*VxI + yI/rI*VyI + zI/rI*VzI;
    CFreal VrI_dimless = xI_dimless/rI_dimless*VxI_dimless + yI_dimless/rI_dimless*VyI_dimless + zI_dimless/rI_dimless*VzI_dimless;
    CFreal VthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VyI_dimless - rhoI_dimless/rI_dimless*VzI_dimless;
    CFreal VphiI_dimless = -yI_dimless/rhoI_dimless*VxI_dimless + xI_dimless/rhoI_dimless*VyI_dimless;
    

  // VrG = 2.*Vr_user - VrI;
  // Extrapolated according to [Pomoell 2012 ApJ]
  //CFreal VrG = VrI;

  // Dirichlet BC for V: Vr = VrParker(r=RSun), Vtheta=0, Vphi=0
 
//  CFreal VrBoundary = 9165.372593288523 - 1.8204851612575266e-05*rBoundary + 1.0434707449784796e-14*pow(rBoundary,2)
//     - 1.4186777294157635e-24*pow(rBoundary,3) + 8.374676041709383e-35*pow(rBoundary,4)
//     - 1.8430711555714095e-45*pow(rBoundary,5);

  CFreal VrBoundary_dimless = (-0.42108+26.794-628.89+6279.7-12899.60846+6026.2007+1445.2)/(2.2e-4/std::sqrt(1.25e-6*1.67e-13));
  CFreal VthetaBoundary_dimless = 0.0; 
  CFreal VphiBoundary_dimless = 0.0; 


    //VrG = densityI*VrI*rI*rI/(densityG*1.67e-13*rG*rG);
    //CFreal VrG_dimless = (*innerState)[0]*VrI_dimless*rI_dimless*rI_dimless/(densityG*rG_dimless*rG_dimless);

    //CFreal VrG_dimless = VrI_dimless;



    //densityG is already dimless. But better give better variable names here!
        //std::cout << "Vr = " << VrG << endl;


    //VthetaG = VrG*BthetaG/BrG; // Make the two vectors parallel to each other 
    // May 22: just a test:
    //VthetaI = xI*zI/(rhoI*rI)*VxI + yI*zI/(rhoI*rI)*VyI - rhoI/rI*VzI;
    //VthetaG = VthetaI;  // extrapol. acc. to Jens
  //VrG = 2.*VrBoundary - VrI;
  CFreal VthetaG_dimless = VthetaI_dimless*(*innerState)[0]/densityG; // <- JensBC //CFreal VthetaG_dimless = - VthetaI_dimless;
  CFreal VphiG_dimless = VphiI_dimless*(*innerState)[0]/densityG; // <- JensBC //CFreal VphiG_dimless = - VphiI_dimless;
  CFreal VrG_dimless = VrI_dimless*(*innerState)[0]/densityG; //<- JensBC //CFreal VrG_dimless = 2*VrBoundary_dimless - VrI_dimless;


  //VxG = xG/rG*VrG - yG/rhoG*VphiG + xG*zG/(rhoG*rG)*VthetaG;
  //VyG = yG/rG*VrG + xG/rhoG*VphiG + yG*zG/(rhoG*rG)*VthetaG;
  //VzG = zG/rG*VrG - rhoG/rG*VthetaG;

  CFreal VxG_dimless = xG_dimless/rG_dimless*VrG_dimless - yG_dimless/rhoG_dimless*VphiG_dimless + xG_dimless*zG_dimless/(rhoG_dimless*rG_dimless)*VthetaG_dimless;
  CFreal VyG_dimless = yG_dimless/rG_dimless*VrG_dimless + xG_dimless/rhoG_dimless*VphiG_dimless + yG_dimless*zG_dimless/(rhoG_dimless*rG_dimless)*VthetaG_dimless;
  CFreal VzG_dimless = zG_dimless/rG_dimless*VrG_dimless - rhoG_dimless/rG_dimless*VthetaG_dimless;

    //VthetaI = VrI*(*innerState)[1]/(*innerState)[0];  // Vr*Btheta/Br
    //VxI = std::sin(thetaI)*std::cos(phiI)*VrI + std::cos(thetaI)*std::cos(phiI)*VthetaI - std::sin(phiI)*VphiI;
    //VyI = std::sin(thetaI)*std::sin(phiI)*VrI + std::cos(thetaI)*std::sin(phiI)*VthetaI + std::cos(phiI)*VphiI;
    //VzI = std::cos(thetaI)*VrI - std::sin(thetaI)*VthetaI;
    // VALUES ARE EXTRAPOLATED TO THE GHOST CONTINUOUSLY (SO NO FIXED DIRICHLET CONDITION):
    //(*ghostState)[1]  = VxG/(2.2e-4/sqrt(1.256e-6*1.67e-13));
    (*ghostState)[1]  = VxG_dimless;
    //std::cout << "ghoststate Vx = " << VxG << endl;
    //(*ghostState)[2] = VyG/(2.2e-4/sqrt(1.256e-6*1.67e-13));
    (*ghostState)[2] = VyG_dimless;
    //(*ghostState)[3] = VzG/(2.2e-4/sqrt(1.256e-6*1.67e-13));
    (*ghostState)[3] = VzG_dimless;
/*
  } else {  // INSIDE THE DEAD-ZONE: DIRICHLET CONDITIONS
      //std::cout << "in the dead-zone" << endl; 

      // Vr and Vtheta are both zero at the boundary
      // VrG = 2*VrB - VrI = -VrI
      // Compute first VrI from VxI, VyI, VzI
      ////////VxI = (*innerState)[9];
      /////////VyI = (*innerState)[10];
      ///////VzI = (*innerState)[11];
      VrI = xI/rI*VxI + yI/rI*VyI + zI/rI*VzI;
      VrG = - VrI;  // i.e. at the equator no inflow


      //VrG = 2*Vr_user- VrI;
      //std::cout << "Vr = " << VrG << endl;
      // May 22: just a test:
      //VrG = densityI*VrI*rI*rI/(densityG*rG*rG);


  VthetaI = xI*zI/(rhoI*rI)*VxI + yI*zI/(rhoI*rI)*VyI - rhoI/rI*VzI;
      VthetaG = -VthetaI;
      // Now compute VxG, VyG, and VzG from VrG, VthetaG, and VphiG:
  VxG = xG/rG*VrG - yG/rhoG*VphiG + xG*zG/(rhoG*rG)*VthetaG;
  VyG = yG/rG*VrG + xG/rhoG*VphiG + yG*zG/(rhoG*rG)*VthetaG;
  VzG = zG/rG*VrG - rhoG/rG*VthetaG;


      //(*ghostState)[8] = 2.*(*_dimState)[8] - (*innerState)[8];

      /////////VrI = 0.0;
      /////////VthetaI = 0.0;
      ////VxI = std::sin(thetaI)*std::cos(phiI)*VrI + std::cos(thetaI)*std::cos(phiI)*VthetaI - std::sin(phiI)*VphiI;
      /////VyI = std::sin(thetaI)*std::sin(phiI)*VrI + std::cos(thetaI)*std::sin(phiI)*VthetaI + std::cos(phiI)*VphiI;
      //////VzI = std::cos(thetaI)*VrI - std::sin(thetaI)*VthetaI;
      (*ghostState)[9]  = VxG;
      //std::cout << "ghoststate Vx = " << (*ghostState)[9] << endl;
      (*ghostState)[10] = VyG;

      (*ghostState)[11] = VzG;
  }
*/










  // The dipole (not the PFSS):
  //CFreal BrB = 2.*B0dip*RSun*RSun*RSun*std::cos(thetaBoundary)/(rBoundary*rBoundary*rBoundary);
  //CFreal BthetaB = B0dip*RSun*RSun*RSun*std::sin(thetaBoundary)/(rBoundary*rBoundary*rBoundary);
/*
  CFreal BthetaG = BthetaI*pow(rI,3)/pow(rG,3);
  CFreal BphiG = BphiI;
  CFreal c = - 1.0*RSun*pow(1./2.5,3)/(2.0 + pow(1./2.5,3));  // B0dip = 2.2e-4 left out -> code units
  CFreal BrPFSSboundary = -2.*c*cos(thetaBoundary)/pow(rBoundary,3)*pow(RSS,2) - c*cos(thetaBoundary)/RSS;
  CFreal BrG = 2.*BrPFSSboundary - BrI;
//(2.*(-2.2e-4*2.5*6.9551e8*pow(1./2.5,3)/(2.+pow(1/2.5,3)))*zBoundary/rBoundary*(-2.*pow(2.5*6.9551e8,2)/pow(rBoundary,3)-1./(2.5*6.9551e8)) - BrI)/2.2e-4;




  CFreal BxB = -xBoundary*zBoundary/pow(rBoundary,2)*(-2.2e-4*2.5*6.9551e8*pow(1./2.5,3)/(2.+pow(1/2.5,3)))*((2.*pow(6.9551e8,2)/pow(rBoundary,3)+1./(2.5*6.9551e8))+1./rBoundary*(pow(2.5*6.9551e8/rBoundary,2)-(rBoundary/(2.5*6.9551e8))));
  CFreal ByB = -yBoundary*zBoundary/pow(rBoundary,2)*(-2.2e-4*2.5*6.9551e8*pow(1./2.5,3)/(2.+pow(1./2.5,3)))*((2.*pow(6.9551e8,2)/pow(rBoundary,3)+1./(2.5*6.9551e8))+1./rBoundary*(pow(2.5*6.9551e8/rBoundary,2)-(rBoundary/(2.5*6.9551e8))));

  CFreal BzB = -pow(zBoundary,2)/pow(rBoundary,2)*(-2.2e-4*2.5*6.9551e8*pow(1./2.5,3)/(2.+pow(1./2.5,3)))*(2.*pow(6.9551e8,2)/pow(rBoundary,3)+1./(2.5*6.9551e8))+(xBoundary*xBoundary+yBoundary*yBoundary)/pow(rBoundary,3)*(-2.2e-4*2.5*6.9551e8*pow(1./2.5,3)/(2.+pow(1./2.5,3)))*(pow(2.5*6.9551e8/rBoundary,2)-(rBoundary/(2.5*6.9551e8)));


  //CFreal BrG = 2.*BrB - BrI;
  //CFreal BthetaG = 2.*BthetaB - BthetaI;
  //CFreal BphiB = 0.; //TRY ALSO THIS AS THE DIPOLE HAS NO PHI COMPONENT AND VAN DER HOLST SETS IT TO ZERO AT THE BOUNDARY, WHILE SKRALAN USES EXTRAPOLATION TO THE GHOST HERE
  //CFreal BphiG = 2.*BphiB - BphiI;
  //CFreal BphiG = BphiI;

  CFreal BxG = xG/rG*BrG - yG/rhoG*BphiG + xG*zG/(rhoG*rG)*BthetaG;
  CFreal ByG = yG/rG*BrG + xG/rhoG*BphiG + yG*zG/(rhoG*rG)*BthetaG;
  CFreal BzG = zG/rG*BrG - rhoG/rG*BthetaG;




  (*ghostState)[0] = 2.*BxB - BxI;
  (*ghostState)[1] = 2.*ByB - ByI;
  (*ghostState)[2] = 2.*BzB - BzI;

  (*ghostState)[4] = BxG;
  (*ghostState)[5] = ByG;
  (*ghostState)[6] = BzG;
*/



  /*CFreal EthetaBoundary = -VphiBoundary*BrBoundary;

  CFreal ErG = - ErI;
  CFreal EthetaG = 2*EthetaBoundary - EthetaI;
  CFreal EphiG = - EphiI;

  CFreal ExG = xG/rG*ErG - yG/rhoG*EphiG + xG*zG/(rhoG*rG)*EthetaG;
  CFreal EyG = yG/rG*ErG + xG/rhoG*EphiG + yG*zG/(rhoG*rG)*EthetaG;
  CFreal EzG = zG/rG*ErG - rhoG/rG*EthetaG;


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

/*
  (*ghostState)[3] = 2.*ExB - (*innerState)[3];
  (*ghostState)[4] = 2.*EyB - (*innerState)[4];
  (*ghostState)[5] = 2.*EzB - (*innerState)[5];
*/



  // Temperature kept constant at 1.5e6 K (set in CFcase)
  // Pressure at inner boundary:
  CFreal PressureBoundary_code_units = (2.0*1.67e-13*1.38e-23*1.5e6/(1.27*1.67e-27))/(pow(2.2e-4,2)/1.2566e-6);

  (*ghostState)[7] = 2.*PressureBoundary_code_units - (*innerState)[7];
  //(*ghostState)[7] = 2.*(*_dimState)[7] - (*innerState)[7];
  //std::cout << "PB_cu = " << PressureBoundary_code_units << endl;
  //std::cout << "PI_cu = " << (*innerState)[7] << endl;
  
 } else {

  // coordinate of the boundary point
  _bCoord = (innerState->getCoordinates() +
             ghostState->getCoordinates());
  _bCoord *= 0.5;

  // ghostState = 2*bcState - innerState
  _vFunction.evaluate(_bCoord, *_dimState);

  if(_inputAdimensionalValues)
  {
    *ghostState = *_dimState;
  }
  else
  {
    _varSet->setAdimensionalValues(*_dimState, *ghostState);
  }
  
  *ghostState *= 2.;
  *ghostState -= *innerState;
  
  //  cf_assert(_projectionIDs.size() > 0);
  for (CFuint i = 0; i < _projectionIDs.size(); ++i) {
    const CFuint varID = _projectionIDs[i];
    (*ghostState)[varID] = (*innerState)[varID];
  }

  // during initialization phase we store the initial solution values to be used as BC 
  if (m_initialSolutionMap.size() > 0) {
    /// map faces to corresponding TRS and index inside that TRS
    SafePtr<MapGeoToTrsAndIdx> mapGeoToTrs =
      MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
    const CFuint faceIdx = mapGeoToTrs->getIdxInTrs(face->getID());
    const string name = getCurrentTRS()->getName();
    SafePtr<RealVector> initialValues = m_initialSolutionMap.find(name);
    const CFuint nbVars = m_initialSolutionIDs.size();
    const CFuint startID = faceIdx*nbVars;
    for (CFuint i = 0; i < nbVars; ++i) {
      const CFuint varID = m_initialSolutionIDs[i];
      const CFuint idx = startID+i;
      cf_assert(idx < initialValues->size());
      (*ghostState)[varID] = 2.*(*initialValues)[idx] - (*innerState)[varID];
      CFLog(DEBUG_MIN, "SuperInletProjection::setGhostState() => [" << varID << "] => " << (*initialValues)[idx] << " | " << (*innerState)[varID] << "\n");
    }
  }

  //CFLog(INFO, "SuperInletProjection::setGhostState() => ghost,dimst,inner = "
  //	<< (*ghostState)[8] << ", " << (*_dimState)[8] << ", e" << (*innerState)[8] << "\n");
  
  /*if ((*ghostState)[8] < 1e-32) {
    CFLog(INFO, "SuperInletProjection::setGhostState() => ghost = " << (*ghostState)[8] <<"\n");
  }

  if ((*innerState)[8] < 1e-32) {
    CFLog(INFO, "SuperInletProjection::setGhostState() => inner = " << (*innerState)[8] <<"\n");
  }

  if ((*_dimState)[8] < 1e-32) {
    CFLog(INFO, "SuperInletProjection::setGhostState() => dimst = " << (*_dimState)[8] <<"\n");
    }*/
}
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
