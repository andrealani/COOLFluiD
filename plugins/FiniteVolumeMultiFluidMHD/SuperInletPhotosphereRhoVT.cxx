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
  

  CFreal xG = 0.;   // Coordinates of the Ghost cell
  CFreal yG = 0.;
  CFreal zG = 0.;
  CFreal rG = 0.;
  CFreal thetaG = 0.;
  CFreal phiG = 0.;
  CFreal latG = 0.;
  CFreal xI = 0.;  // Coordinates of the inner cell
  CFreal yI = 0.;
  CFreal zI = 0.;
  CFreal rI = 0.;
  CFreal thetaI = 0.;
  CFreal phiI = 0.;
  //CFreal latI = 0.;   // will probably be needed when implementing differential rotation
  const CFreal PI = 3.14159265358979;

  xG = ghostState->getCoordinates()[XX];
  yG = ghostState->getCoordinates()[YY];
  zG = ghostState->getCoordinates()[ZZ];
  xI = innerState->getCoordinates()[XX];
  yI = innerState->getCoordinates()[YY];
  zI = innerState->getCoordinates()[ZZ];


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
  rG = std::sqrt(xG*xG + yG*yG + zG*zG);
  //std::cout << "rG = " << rG << endl;
  thetaG = std::atan2(std::sqrt(xG*xG + yG*yG),zG);
  //std::cout << "thetaG = " << thetaG << endl;
  phiG = std::atan2(yG,xG);



  //std::cout << "xI = " << xI << endl;
  //std::cout << "yI = " << yI << endl;
  //std::cout << "zI = " << zI << endl;
  rI = std::sqrt(xI*xI + yI*yI + zI*zI);
  //std::cout << "rI = " << rI << endl;
  thetaI = std::atan2(std::sqrt(xI*xI + yI*yI),zI);
  //std::cout << "thetaI = " << thetaI << endl;
  phiI = std::atan2(yI,xI);
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
 //(*ghostState)[3] = 2.*((*_dimState)[11]*(*_dimState)[1] - (*_dimState)[2]*(*_dimState)[10]) - ((*innerState)[11]*(*innerState)[1] - (*innerState)[2]*(*innerState)[10]); //Ex
  //(*ghostState)[4] = 2.*((*_dimState)[9]*(*_dimState)[2] - (*_dimState)[11]*(*_dimState)[0]) - ((*innerState)[9]*(*innerState)[2] - (*innerState)[11]*(*innerState)[0]); //Ex
  //(*ghostState)[5] = 2.*((*_dimState)[10]*(*_dimState)[0] - (*_dimState)[9]*(*_dimState)[1]) - ((*innerState)[10]*(*innerState)[0] - (*innerState)[9]*(*innerState)[1]); //Ex
  // ... the last 3 lines lead immediately to a PETSC error
  (*ghostState)[3] = 2.*((*_dimState)[11]*(*_dimState)[1] - (*_dimState)[2]*(*_dimState)[10]) - (*innerState)[3]; //Ex
  (*ghostState)[4] = 2.*((*_dimState)[9]*(*_dimState)[2] - (*_dimState)[11]*(*_dimState)[0]) - (*innerState)[4]; //Ey
  (*ghostState)[5] = 2.*((*_dimState)[10]*(*_dimState)[0] - (*_dimState)[9]*(*_dimState)[1]) - (*innerState)[5]; //Ex
// Perfectly conducting wall BC
  //(*ghostState)[3] = - (*innerState)[3] + 2.0*((*innerState)[3]*std::sin(thetaI)*std::cos(phiI) + (*innerState)[4]*std::sin(thetaI)*std::sin(phiI) + (*innerState)[5]*std::cos(thetaI))*std::sin(thetaI)*std::cos(phiI);
  //(*ghostState)[4] = - (*innerState)[4] + 2.0*((*innerState)[3]*std::sin(thetaI)*std::cos(phiI) + (*innerState)[4]*std::sin(thetaI)*std::sin(phiI) + (*innerState)[5]*std::cos(thetaI))*std::sin(thetaI)*std::sin(phiI);
  //(*ghostState)[5] = - (*innerState)[5] + 2.0*((*innerState)[3]*std::sin(thetaI)*std::cos(phiI) + (*innerState)[4]*std::sin(thetaI)*std::sin(phiI) + (*innerState)[5]*std::cos(thetaI))*std::cos(thetaI);
  //(*ghostState)[3] = 0.;
  //(*ghostState)[4] = 0.;
  //(*ghostState)[5] = 0.;

  // Continuous condition for Phi, Psi
  //(*ghostState)[6] = (*innerState)[6];
  //(*ghostState)[7] = (*innerState)[7];

  (*ghostState)[6] = (*innerState)[6];
  (*ghostState)[7] = - (*innerState)[7];


  (*ghostState)[8] = 2.*(*_dimState)[8] - (*innerState)[8];
  //std::cout << "ghoststate rho = " << (*ghostState)[8] << endl;
  
  // VELOCITY: This includes differential rotation but the first version should not
  //CFreal OmegaG = 0.;
  //CFreal OmegaI = 0.;
  //CFreal A = 14.713;   // values from Wikipedia; find a better and more accurate source...
  //CFreal B = -2.396;
  //CFreal C = -1.787;
  
  //OmegaG = A + B*std::sin(latG)*std::sin(latG) + C*std::sin(latG)*std::sin(latG)*std::sin(latG)*std::sin(latG);
  //OmegaI = A + B*std::sin(latI)*std::sin(latI) + C*std::sin(latI)*std::sin(latI)*std::sin(latI)*std::sin(latI);
  // Convert OmegaG and OmegaI from deg/day to rad/s
  // deg/day = (PI/180.*rad)/(24*60*60*sec)
  //OmegaG = OmegaG*(PI/180)/(24*60*60);
  //OmegaI = OmegaI*(PI/180)/(24*60*60);


  CFreal VrG = 0.;
  CFreal VrI = 0.;
  CFreal VthetaG = 0.;
  CFreal VthetaI = 0.;
  CFreal VphiG = 0.;
  CFreal VphiI = 0.;
  CFreal VxG = 0.;
  CFreal VyG = 0.;
  CFreal VzG = 0.;
  CFreal VxI = 0.;
  CFreal VyI = 0.;
  CFreal VzI = 0.;
    
  CFreal BxG = (*ghostState)[0];
  //cout << "Bxfield = " << BxG << endl;
  CFreal ByG = (*ghostState)[1];
  //cout << "Byfield = " << BxG << endl;
  CFreal  BzG = (*ghostState)[2];
  //cout << "Bzfield = " << BxG << endl;
  
  // Convert the B-field at the ghost state into cartesian vector components:
  CFreal BrG = std::sin(thetaG)*std::cos(phiG)*BxG + std::sin(thetaG)*std::sin(phiG)*ByG + std::cos(thetaG)*BzG;
  CFreal BthetaG = std::cos(thetaG)*std::cos(phiG)*BxG + std::cos(thetaG)*std::sin(phiG)*ByG - std::sin(thetaG)*BzG;
  //CFreal BphiG = -std::sin(phiG)*BxG + std::cos(phiG)*ByG;
  
  // VphiG = 2*VphiB - VphiI = -VphiI; // only preliminary before taking differential rotation into account; this value is independent of the zone
  VxI = (*innerState)[9];
  VyI = (*innerState)[10];
  VzI = (*innerState)[11];
  VphiI = -std::sin(phiI)*VxI + std::cos(phiI)*VyI;
  VphiG = -VphiI;

  if (latG  > 0.3927) {  // this value corresponds to 45/2 deg in rad
    //cout << "in the pole region" << endl;
    VrI = std::sin(thetaI)*std::cos(phiI)*(*innerState)[9] + std::sin(thetaI)*std::sin(phiI)*(*innerState)[10] + std::cos(thetaI)*(*innerState)[11];
    VrG = (*innerState)[8]*VrI*rI*rI/((*ghostState)[8]*rG*rG);
    
    
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
    (*ghostState)[9]  = VxG;
    //std::cout << "ghoststate Vx = " << VxG << endl;
    (*ghostState)[10] = VyG;
    (*ghostState)[11] = VzG;
  } else {  // INSIDE THE DEAD-ZONE: DIRICHLET CONDITIONS
      //std::cout << "in the dead-zone" << endl; 

      // Vr and Vtheta are both zero at the boundary
      // VrG = 2*VrB - VrI = -VrI
      // Compute first VrI from VxI, VyI, VzI
      VxI = (*innerState)[9];
      VyI = (*innerState)[10];
      VzI = (*innerState)[11];
      VrI = std::sin(thetaI)*std::cos(phiI)*VxI + std::sin(thetaI)*std::sin(phiI)*VyI + std::cos(thetaI)*VzI;
      


      //VrG = - VrI;
      // May 22: just a test:
      VrG = (*innerState)[8]*VrI*rI*rI/((*ghostState)[8]*rG*rG);


      VthetaI = std::cos(thetaI)*std::cos(phiI)*VxI + std::cos(thetaI)*std::sin(phiI)*VyI - std::sin(thetaI)*VzI;
      VthetaG = -VthetaI;
      // Now compute VxG, VyG, and VzG from VrG, VthetaG, and VphiG:
      VxG = std::sin(thetaG)*std::cos(phiG)*VrG + std::cos(thetaG)*std::cos(phiG)*VthetaG - std::sin(phiG)*VphiG;
      VyG = std::sin(thetaG)*std::sin(phiG)*VrG + std::cos(thetaG)*std::sin(phiG)*VthetaG + std::cos(phiG)*VphiG;
      VzG = std::cos(thetaG)*VrG - std::sin(thetaG)*VthetaG;

      //(*ghostState)[8] = 2.*(*_dimState)[8] - (*innerState)[8];

      VrI = 0.0;
      VthetaI = 0.0;
      VxI = std::sin(thetaI)*std::cos(phiI)*VrI + std::cos(thetaI)*std::cos(phiI)*VthetaI - std::sin(phiI)*VphiI;
      VyI = std::sin(thetaI)*std::sin(phiI)*VrI + std::cos(thetaI)*std::sin(phiI)*VthetaI + std::cos(phiI)*VphiI;
      VzI = std::cos(thetaI)*VrI - std::sin(thetaI)*VthetaI;
      (*ghostState)[9]  = 2.*(*_dimState)[9] - VxI;   // dimState[9] set in CFcase
      //std::cout << "ghoststate Vx = " << (*ghostState)[9] << endl;
      (*ghostState)[10] = 2.*(*_dimState)[10] - VyI;   // dimState[9] set in CFcase
      (*ghostState)[11] = 2.*(*_dimState)[9] - VzI;   // dimState[9] set in CFcase
  }

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
