#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SuperInletRhoVT.hh"
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

MethodCommandProvider<SuperInletRhoVT, CellCenterFVMData,
		      FiniteVolumeNavierStokesModule>
superInletRhoVTFVMCCProvider("SuperInletRhoVTFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SuperInletRhoVT::defineConfigOptions(Config::OptionList& options)
{
  //options.addConfigOption< std::vector<CFreal> >("Ttot","total temperature");
}

//////////////////////////////////////////////////////////////////////////////

SuperInletRhoVT::SuperInletRhoVT(const std::string& name) :
  SuperInletProjection(name)
{
   addConfigOptionsTo(this);
   //_tTotal = std::vector<CFreal>();
   //setParameter("Ttot",&_tTotal);
}

//////////////////////////////////////////////////////////////////////////////

SuperInletRhoVT::~SuperInletRhoVT()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletRhoVT::setGhostState(GeometricEntity *const face)
{
  // first call the base setGhostState()
  SuperInletProjection::setGhostState(face);
  
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  // then overwrite the ghost variables that you need to overwrite
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  CFreal latG = 0.;
  //CFreal latI = 0.;   // will probably be needed when implementing differential rotation
  const CFreal PI = MathTools::MathConsts::CFrealPi();
  const CFreal xG = ghostState->getCoordinates()[XX];
  const CFreal yG = ghostState->getCoordinates()[YY];
  const CFreal zG = ghostState->getCoordinates()[ZZ];
  const CFreal xI = innerState->getCoordinates()[XX];
  const CFreal yI = innerState->getCoordinates()[YY];
  const CFreal zI = innerState->getCoordinates()[ZZ];
  CFreal xBoundary = (xG + xI)/2.0;
  CFreal yBoundary = (yG + yI)/2.0;
  CFreal zBoundary = (zG + zI)/2.0;
  CFreal rBoundary = std::sqrt(xBoundary*xBoundary + yBoundary*yBoundary + zBoundary*zBoundary);
  CFreal thetaBoundary = std::atan2(std::sqrt(xBoundary*xBoundary + yBoundary*yBoundary),zBoundary);
  CFreal phiBoundary = std::atan2(yBoundary,xBoundary);
  
  // Convert boundary values provided in the CFcase file into spherical coordinates
  CFreal Vx_user = (*_dimState)[1];
  CFreal Vy_user = (*_dimState)[2];
  CFreal Vz_user = (*_dimState)[3];
  CFreal Vr_user = std::sin(thetaBoundary)*std::cos(phiBoundary)*Vx_user + std::sin(thetaBoundary)*std::sin(phiBoundary)*Vy_user + std::cos(thetaBoundary)*Vz_user;
  CFreal Vtheta_user = std::cos(thetaBoundary)*std::cos(phiBoundary)*Vx_user + std::cos(thetaBoundary)*std::sin(phiBoundary)*Vy_user - std::sin(thetaBoundary)*Vz_user;
  CFreal Vphi_user = -std::sin(phiBoundary)*Vx_user + std::cos(phiBoundary)*Vy_user;
  
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
  
  //std::cout << "xG = " << xG << endl;
  //std::cout << "yG = " << yG << endl;
  //std::cout << "zG = " << zG << endl;
  const CFreal rG = std::sqrt(xG*xG + yG*yG + zG*zG);
  //std::cout << "rG = " << rG << endl;
  const CFreal thetaG = std::atan2(std::sqrt(xG*xG + yG*yG),zG);
  //std::cout << "thetaG = " << thetaG << endl;
  const CFreal phiG = std::atan2(yG,xG);
  
  //std::cout << "xI = " << xI << endl;
  //std::cout << "yI = " << yI << endl;
  //std::cout << "zI = " << zI << endl;
  const CFreal rI = std::sqrt(xI*xI + yI*yI + zI*zI);
  //std::cout << "rI = " << rI << endl;
  const CFreal thetaI = std::atan2(std::sqrt(xI*xI + yI*yI),zI);
  //std::cout << "thetaI = " << thetaI << endl;
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
//  if (thetaI > -PI && thetaI < -PI*0.5) {
//     latI = std::abs(thetaI) - PI*0.5;
//  } else if (thetaI > -PI*0.5 && thetaI < 0) {
//     latI = PI*0.5 - std::abs(thetaI);
//  } else if (thetaG > 0. && thetaG < PI*0.5) {
//     latI = PI*0.5 - thetaI;
//  } else if (thetaI > PI*0.5 && thetaI < PI) {
//     latI = thetaI - PI*0.5;
//  } else {
//    std::cout << "Error: value of theta for the point in question outside expected range" << endl;
//  }
  
  // DENSITY kept constant (Dirichlet condition)
  
  (*ghostState)[0] = 2.*(*_dimState)[0] - (*innerState)[0];
  
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
  CFreal rhoG = (*ghostState)[0];
  
  // VphiG = 2*VphiB - VphiI = -VphiI; // only preliminary before taking differential rotation into account; this value is independent of the zone
  CFreal rhoI = (*innerState)[0];
  CFreal VxI = (*innerState)[1];
  CFreal VyI = (*innerState)[2];
  CFreal VzI = (*innerState)[3];
  CFreal VphiI = -std::sin(phiI)*VxI + std::cos(phiI)*VyI;
  CFreal VphiG = - VphiI; //2*Omega*rBoundary*std::sin(thetaBoundary) - VphiI;
  //std::cout << "Vphi = " << VphiG << endl;
  
  if (latG  > 0.3927) {  // this value corresponds to 45/2 deg in rad
    //cout << "in the pole region" << endl;
    VrI = std::sin(thetaI)*std::cos(phiI)*VxI + std::sin(thetaI)*std::sin(phiI)*VzI + std::cos(thetaI)*VzI;
    VrG = rhoI*VrI*rI*rI/(rhoG*rG*rG);
    //std::cout << "Vr = " << VrG << endl;
    
    
    //VthetaG = VrG*BthetaG/BrG; // Make the two vectors parallel to each other 
    // May 22: just a test:
    VthetaI = std::cos(thetaI)*std::cos(phiI)*VxI + std::cos(thetaI)*std::sin(phiI)*VyI - std::sin(thetaI)*VzI;
    VthetaG = -VthetaI;
    
    
    VxG = std::sin(thetaG)*std::cos(phiG)*VrG + std::cos(thetaG)*std::cos(phiG)*VthetaG - std::sin(phiG)*VphiG;
    VyG = std::sin(thetaG)*std::sin(phiG)*VrG + std::cos(thetaG)*std::sin(phiG)*VthetaG + std::cos(phiG)*VphiG;
    VzG = std::cos(thetaG)*VrG - std::sin(thetaG)*VthetaG;
    //VthetaI = VrI*(*innerState)[1]/(*innerState)[0];  // Vr*Btheta/Br
    //VxI = std::sin(thetaI)*std::cos(phiI)*VrI + std::cos(thetaI)*std::cos(phiI)*VthetaI - std::sin(phiI)*VphiI;
    //VyI = std::sin(thetaI)*std::sin(phiI)*VrI + std::cos(thetaI)*std::sin(phiI)*VthetaI + std::cos(phiI)*VphiI;
    //VzI = std::cos(thetaI)*VrI - std::sin(thetaI)*VthetaI;
    // VALUES ARE EXTRAPOLATED TO THE GHOST CONTINUOUSLY (SO NO FIXED DIRICHLET CONDITION):
    (*ghostState)[1]  = VxG;
    //std::cout << "ghoststate Vx = " << VxG << endl;
    (*ghostState)[2] = VyG;
    (*ghostState)[3] = VzG;
  } else {  // INSIDE THE DEAD-ZONE: DIRICHLET CONDITIONS
    //std::cout << "in the dead-zone" << endl; 
    
    // Vr and Vtheta are both zero at the boundary
    // VrG = 2*VrB - VrI = -VrI
    // Compute first VrI from VxI, VyI, VzI
    ////////VxI = (*innerState)[9];
    /////////VyI = (*innerState)[10];
    ///////VzI = (*innerState)[11];
    VrI = std::sin(thetaI)*std::cos(phiI)*VxI + std::sin(thetaI)*std::sin(phiI)*VyI + std::cos(thetaI)*VzI;
    
    
    
    //VrG = 2*Vr_user- VrI;
    //std::cout << "Vr = " << VrG << endl;
    // May 22: just a test:
    VrG = rhoI*VrI*rI*rI/(rhoG*rG*rG);
    
    
    VthetaI = std::cos(thetaI)*std::cos(phiI)*VxI + std::cos(thetaI)*std::sin(phiI)*VyI - std::sin(thetaI)*VzI;
    VthetaG = -VthetaI;
    // Now compute VxG, VyG, and VzG from VrG, VthetaG, and VphiG:
    VxG = std::sin(thetaG)*std::cos(phiG)*VrG + std::cos(thetaG)*std::cos(phiG)*VthetaG - std::sin(phiG)*VphiG;
    VyG = std::sin(thetaG)*std::sin(phiG)*VrG + std::cos(thetaG)*std::sin(phiG)*VthetaG + std::cos(phiG)*VphiG;
    VzG = std::cos(thetaG)*VrG - std::sin(thetaG)*VthetaG;
    
    //(*ghostState)[8] = 2.*(*_dimState)[8] - (*innerState)[8];
    
    /////////VrI = 0.0;
    /////////VthetaI = 0.0;
    ////VxI = std::sin(thetaI)*std::cos(phiI)*VrI + std::cos(thetaI)*std::cos(phiI)*VthetaI - std::sin(phiI)*VphiI;
    /////VyI = std::sin(thetaI)*std::sin(phiI)*VrI + std::cos(thetaI)*std::sin(phiI)*VthetaI + std::cos(phiI)*VphiI;
    //////VzI = std::cos(thetaI)*VrI - std::sin(thetaI)*VthetaI;
    (*ghostState)[1]  = VxG;
    //std::cout << "ghoststate Vx = " << (*ghostState)[9] << endl;
    (*ghostState)[2] = VyG;
    (*ghostState)[3] = VzG;
  }
  
  // Temperature kept constant at 1.5e6 K (set in CFcase)
  (*ghostState)[4] = 2.*(*_dimState)[4] - (*innerState)[4];
}
      
//////////////////////////////////////////////////////////////////////////////

void SuperInletRhoVT::setup()
{
  SuperInletProjection::setup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
