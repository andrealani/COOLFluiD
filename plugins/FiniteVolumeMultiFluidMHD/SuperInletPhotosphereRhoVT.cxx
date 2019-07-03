#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/SuperInletPhotosphereRhoVT.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalConsts.hh"

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

MethodCommandProvider<SuperInletPhotosphereRhoVT, CellCenterFVMData,
		      FiniteVolumeMultiFluidMHDModule>
superInletPhotosphereRhoVTFVMCCProvider("SuperInletPhotosphereRhoVTFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SuperInletPhotosphereRhoVT::defineConfigOptions(Config::OptionList& options)
{
  //options.addConfigOption< std::vector<CFreal> >("Ttot","total temperature");
}

//////////////////////////////////////////////////////////////////////////////

SuperInletPhotosphereRhoVT::SuperInletPhotosphereRhoVT(const std::string& name) :
  SuperInletProjection(name)
{
   addConfigOptionsTo(this);
   //_tTotal = std::vector<CFreal>();
   //setParameter("Ttot",&_tTotal);
}

//////////////////////////////////////////////////////////////////////////////

SuperInletPhotosphereRhoVT::~SuperInletPhotosphereRhoVT()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletPhotosphereRhoVT::setGhostState(GeometricEntity *const face)
{
  // first call the base setGhostState()
  SuperInletProjection::setGhostState(face);

  if (m_initialSolutionMap.size()>0) {
  
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // then overwrite the ghost variables that you need to overwrite
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  //cout << "----------------------------------------------------" << endl;
  CFreal latG = 0.;
  //CFreal latI = 0.;   // will probably be needed when implementing differential rotation
  const CFreal PI = MathTools::MathConsts::CFrealPi();
  const CFreal xG = ghostState->getCoordinates()[XX];
  //cout << "xG = " << xG << endl;
  const CFreal yG = ghostState->getCoordinates()[YY];
  //cout << "yG = " << yG << endl;
  const CFreal zG = ghostState->getCoordinates()[ZZ];
  //cout << "zG = " << zG << endl;
  const CFreal xI = innerState->getCoordinates()[XX];
 //cout << "xI = " << xI << endl;
  const CFreal yI = innerState->getCoordinates()[YY];
 //cout << "yI = " << yI << endl;
  const CFreal zI = innerState->getCoordinates()[ZZ];
 //cout << "zI = " << zI << endl;
  CFreal xBoundary = (xG + xI)/2.0;
 //cout << "xB = " << xBoundary << endl;
  CFreal yBoundary = (yG + yI)/2.0;
 //cout << "yB = " << yBoundary << endl;
  CFreal zBoundary = (zG + zI)/2.0;
 //cout << "zB = " << zBoundary << endl;
  CFreal rBoundary = std::sqrt(xBoundary*xBoundary + yBoundary*yBoundary + zBoundary*zBoundary);
 //cout << "rB = " << rBoundary << endl;
  CFreal rhoBoundary = std::sqrt(xBoundary*xBoundary + yBoundary*yBoundary);
  CFreal thetaBoundary = std::atan2(std::sqrt(xBoundary*xBoundary + yBoundary*yBoundary),zBoundary);
 //cout << "thetaB = " << thetaBoundary << endl;
  CFreal phiBoundary = std::atan2(yBoundary,xBoundary);
 //cout << "phiB = " << phiBoundary << endl;



  // Convert boundary values provided in the CFcase file into spherical coordinates
  CFreal Vx_user = (*_dimState)[9];
  CFreal Vy_user = (*_dimState)[10];
  CFreal Vz_user = (*_dimState)[11];
  CFreal Vr_user = xBoundary/rBoundary*Vx_user + yBoundary/rBoundary*Vy_user + zBoundary/rBoundary*Vz_user;
  CFreal Vtheta_user = xBoundary*zBoundary/(rhoBoundary*rBoundary)*Vx_user + yBoundary*zBoundary/(rhoBoundary*rBoundary)*Vy_user - rhoBoundary/rBoundary*Vz_user;
  CFreal Vphi_user = -yBoundary/rhoBoundary*Vx_user + xBoundary/rhoBoundary*Vy_user;





  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = 0.;
  CFreal sum2 = nx*nx + ny*ny;
  if (dim == DIM_3D) {
    nz = normals[startID + 2];
    sum2 += nz*nz;
  }
  const CFreal invFaceLength = 1./sqrt(sum2);
  nx *= invFaceLength;
  ny *= invFaceLength;
  nz *= invFaceLength;
  
  // here the magic happen... Examples:
  // (*ghostState)[varN] = (*innerState)[varN]; // Neumann
  // (*ghostState)[varD] = 2.*inlet_varD - (*innerState)[varD]; // Dirichlet
  // where inlet_varD = (*_dimState)[varN]
  // Set all flow variables here rho = 8, Vx, Vy, Vz (9, 10, 11), T (12)

  const CFreal rG = std::sqrt(xG*xG + yG*yG + zG*zG);
 //cout << "rG = " << rG << endl;
  const CFreal rhoG = std::sqrt(xG*xG + yG*yG);

  const CFreal thetaG = std::atan2(std::sqrt(xG*xG + yG*yG),zG);
 //cout << "thetaG = " << thetaG << endl;
  const CFreal phiG = std::atan2(yG,xG);
 //cout << "phiG = " << phiG << endl;



  const CFreal rI = std::sqrt(xI*xI + yI*yI + zI*zI);
  const CFreal thetaI = std::atan2(std::sqrt(xI*xI + yI*yI),zI);
  const CFreal rhoI = std::sqrt(xI*xI + yI*yI);
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
    std::cout << "thetaG = " << thetaG << endl;

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

  (*ghostState)[6] = -(*innerState)[6]; // a + in Alejandro's case
  (*ghostState)[7] = -(*innerState)[7]; // - in Alejandro's case

  CFreal densityG = 2.*(*_dimState)[8] - (*innerState)[8];
  (*ghostState)[8] = densityG;
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


  CFreal VrG = 0.;
  CFreal VrI = 0.;
  CFreal VthetaG = 0.;
  CFreal VthetaI = 0.;
  CFreal VxG = 0.;
  CFreal VyG = 0.;
  CFreal VzG = 0.;




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

  CFreal BxI = (*innerState)[0];
 //cout << "BxI = " << BxI << endl;
  CFreal ByI = (*innerState)[1];
 //cout << "ByI = " << ByI << endl;
  CFreal BzI = (*innerState)[2];
 //cout << "BzI = " << BzI << endl;
  CFreal BrI = xI/rI*BxI + yI/rI*ByI + zI/rI*BzI;
 //cout << "BrI = " << BrI << endl;
  CFreal BthetaI = xI*zI/(rhoI*rI)*BxI + yI*zI/(rhoI*rI)*ByI - rhoI/rI*BzI;
 //cout << "BthetaI = " << BthetaI << endl;
  CFreal BphiI = -yI/rhoI*BxI + xI/rhoI*ByI;
 //cout << "BphiI = " << BphiI << endl;



  // VphiG = 2*VphiB - VphiI = -VphiI; // only preliminary before taking differential rotation into account; this value is independent of the zone
  CFreal densityI = (*innerState)[8];
  CFreal VxI = (*innerState)[9];
  CFreal VyI = (*innerState)[10];
  CFreal VzI = (*innerState)[11];
  CFreal VphiI = -std::sin(phiI)*VxI + std::cos(phiI)*VyI;
  CFreal VphiBoundary = Omega*rBoundary*std::sin(thetaBoundary);
  //CFreal VphiG = 2*VphiBoundary - VphiI; //- VphiI; //
  CFreal VphiG = -VphiI;  

//std::cout << "Vphi = " << VphiG << endl;

//  if (latG  > 0.3927) {  // this value corresponds to 45/2 deg in rad
    //cout << "in the pole region" << endl;
    VrI = xI/rI*VxI + yI/rI*VyI + zI/rI*VzI;
    VrG = 2.*Vr_user - VrI;
    //VrG = densityI*VrI*rI*rI/(densityG*rG*rG);
        //std::cout << "Vr = " << VrG << endl;


    //VthetaG = VrG*BthetaG/BrG; // Make the two vectors parallel to each other 
    // May 22: just a test:
    VthetaI = xI*zI/(rhoI*rI)*VxI + yI*zI/(rhoI*rI)*VyI - rhoI/rI*VzI;
    VthetaG = -VthetaI;


  VxG = xG/rG*VrG - yG/rhoG*VphiG + xG*zG/(rhoG*rG)*VthetaG;
  VyG = yG/rG*VrG + xG/rhoG*VphiG + yG*zG/(rhoG*rG)*VthetaG;
  VzG = zG/rG*VrG - rhoG/rG*VthetaG;
    //VthetaI = VrI*(*innerState)[1]/(*innerState)[0];  // Vr*Btheta/Br
    //VxI = std::sin(thetaI)*std::cos(phiI)*VrI + std::cos(thetaI)*std::cos(phiI)*VthetaI - std::sin(phiI)*VphiI;
    //VyI = std::sin(thetaI)*std::sin(phiI)*VrI + std::cos(thetaI)*std::sin(phiI)*VthetaI + std::cos(phiI)*VphiI;
    //VzI = std::cos(thetaI)*VrI - std::sin(thetaI)*VthetaI;
    // VALUES ARE EXTRAPOLATED TO THE GHOST CONTINUOUSLY (SO NO FIXED DIRICHLET CONDITION):
    (*ghostState)[9]  = VxG;
    //std::cout << "ghoststate Vx = " << VxG << endl;
    (*ghostState)[10] = VyG;
    (*ghostState)[11] = VzG;
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


  ////CFreal ExG = (*ghostState)[3];
  ////CFreal EyG = (*ghostState)[4];
  ////CFreal EzG = (*ghostState)[5];
  CFreal ExI = (*innerState)[3];
  CFreal EyI = (*innerState)[4];
  CFreal EzI = (*innerState)[5];

  /////CFreal ErG = xG/rG*ExG + yG/rG*EyG + zG/rG*EzG;
  /////CFreal EthetaG = xG*zG/(rhoG*rG)*ExG + yG*zG/(rhoG*rG)*EyG - rhoG/rG*EzG;
  ////CFreal EphiG = -yG/rhoG*ExG + xG/rhoG*EyG;
  CFreal ErI = xI/rI*ExI + yI/rI*EyI + zI/rI*EzI;
  CFreal EthetaI = xI*zI/(rhoI*rI)*ExI + yI*zI/(rhoI*rI)*EyI - rhoI/rI*EzI;
  CFreal EphiI = -yI/rhoI*ExI + xI/rhoI*EyI;

//  CFreal BxBoundary = (BxG + BxI)/2.;
// //cout << "BxB  " << BxBoundary << endl;
//  CFreal ByBoundary = (ByG + ByI)/2.;
// //cout << "ByB = " << ByBoundary << endl;
//  CFreal BzBoundary = (BzG + BzI)/2.;
// //cout << "BzB = " << BzBoundary << endl;


  CFreal RSun = 6.9551e8; // m
  CFreal B0dip = 0.00022; // Tesla
  CFreal BrB = 2.*B0dip*RSun*RSun*RSun*cos(thetaBoundary)/(rBoundary*rBoundary*rBoundary);
  CFreal BthetaB = B0dip*RSun*RSun*RSun*sin(thetaBoundary)/(rBoundary*rBoundary*rBoundary);
  CFreal BrG = 2.*BrB - BrI;
  CFreal BthetaG = 2.*BthetaB - BthetaI;
  //CFreal BphiB = 0.; //TRY ALSO THIS AS THE DIPOLE HAS NO PHI COMPONENT AND VAN DER HOLST SETS IT TO ZERO AT THE BOUNDARY, WHILE SKRALAN USES EXTRAPOLATION TO THE GHOST HERE
  //CFreal BphiG = 2.*BphiB - BphiI;
  CFreal BphiG = BphiI;

  CFreal BxG = xG/rG*BrG - yG/rhoG*BphiG + xG*zG/(rhoG*rG)*BthetaG;
  CFreal ByG = yG/rG*BrG + xG/rhoG*BphiG + yG*zG/(rhoG*rG)*BthetaG;
  CFreal BzG = zG/rG*BrG - rhoG/rG*BthetaG;

  (*ghostState)[0] = BxG; //Ex
  (*ghostState)[1] = ByG; //Ey
  (*ghostState)[2] = BzG;

  /*CFreal EthetaBoundary = -VphiBoundary*BrBoundary;

  CFreal ErG = - ErI;
  CFreal EthetaG = 2*EthetaBoundary - EthetaI;
  CFreal EphiG = - EphiI;

  CFreal ExG = xG/rG*ErG - yG/rhoG*EphiG + xG*zG/(rhoG*rG)*EthetaG;
  CFreal EyG = yG/rG*ErG + xG/rhoG*EphiG + yG*zG/(rhoG*rG)*EthetaG;
  CFreal EzG = zG/rG*ErG - rhoG/rG*EthetaG;

  */

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

  (*ghostState)[3] = 2.*ExB - (*innerState)[3];
  (*ghostState)[4] = 2.*EyB - (*innerState)[4];
  (*ghostState)[5] = 2.*EzB - (*innerState)[5];



  // Temperature kept constant at 1.5e6 K (set in CFcase)
  (*ghostState)[12] = 2.*(*_dimState)[12] - (*innerState)[12];
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void SuperInletPhotosphereRhoVT::setup()
{
  SuperInletProjection::setup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
