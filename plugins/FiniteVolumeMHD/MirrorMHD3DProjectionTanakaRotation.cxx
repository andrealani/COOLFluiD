#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MirrorMHD3DProjectionTanakaRotation.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "Framework/MeshData.hh"
  
//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

MethodCommandProvider<MirrorMHD3DProjectionTanakaRotation, CellCenterFVMData, FiniteVolumeMHDModule> 
mirrorMHD3DProjectionTanakaRotationFVMCCProvider("MirrorMHD3DProjectionTanakaRotationFVMCC");

//////////////////////////////////////////////////////////////////////
   
void MirrorMHD3DProjectionTanakaRotation::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("rhoFixed","Density value that is to be fixed when the normal velocity is outward from the body.");
   options.addConfigOption< CFreal >("T","Non-dimensional period of rotation of the planet/satellite.");
}

//////////////////////////////////////////////////////////////////////

MirrorMHD3DProjectionTanakaRotation::MirrorMHD3DProjectionTanakaRotation(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
  addConfigOptionsTo(this);

  _rhoFixed = 1.0;
  setParameter("rhoFixed",&_rhoFixed);

  _T = 1.0;
  setParameter("T",&_T);
}

//////////////////////////////////////////////////////////////////////

MirrorMHD3DProjectionTanakaRotation::~MirrorMHD3DProjectionTanakaRotation()
{
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionTanakaRotation::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionTanakaRotation::setup()
{
  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionTanakaRotation::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
  
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = normals[startID + 2];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx *= invFaceLength;
  ny *= invFaceLength;
  nz *= invFaceLength;

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  /////////////////////////////////////////////////////////////////////////////////

  // boundary conditions for rho, VX, VY and VZ according to Tanaka's JCP, Vol.111,
  // pp.381-389, 1994 paper, which were implemented prior to Gabor Toth's suggestions

 //  const CFreal vn = _dataInnerState[MHDTerm::VX]*nx + 
//     _dataInnerState[MHDTerm::VY]*ny + _dataInnerState[MHDTerm::VZ]*nz;
//   const CFreal bn = _dataInnerState[MHDTerm::BX]*nx +
//     _dataInnerState[MHDTerm::BY]*ny + _dataInnerState[MHDTerm::BZ]*nz;

//   if (vn > 0.0)
//     _dataGhostState[MHDTerm::RHO] = _rhoFixed;
//   else
//     _dataGhostState[MHDTerm::RHO] = _dataInnerState[MHDTerm::RHO];

//   RealVector BDipoleInner(PhysicalModelStack::getActive()->getDim()), 
//     BDipoleGhost(PhysicalModelStack::getActive()->getDim());

//   _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);

//   BDipoleInner = _varSet->getMagneticDipole(*innerState);
//   BDipoleGhost = _varSet->getMagneticDipole(*ghostState);
//   const CFreal Bx0 = 0.5*(BDipoleInner[0]+BDipoleGhost[0]);
//   const CFreal By0 = 0.5*(BDipoleInner[1]+BDipoleGhost[1]);
//   const CFreal Bz0 = 0.5*(BDipoleInner[2]+BDipoleGhost[2]);
//   const CFreal mx = _dataInnerState[MHDTerm::RHO]*_dataInnerState[MHDTerm::VX];
//   const CFreal my = _dataInnerState[MHDTerm::RHO]*_dataInnerState[MHDTerm::VY];
//   const CFreal mz = _dataInnerState[MHDTerm::RHO]*_dataInnerState[MHDTerm::VZ];    
//   const CFreal mdotB0 = mx*Bx0 + my*By0 + mz*Bz0;
//   const CFreal alpha = mdotB0/(Bx0*Bx0+By0*By0+Bz0*Bz0);
  
//   const CFreal mxGhost = 2.0*alpha*Bx0 - mx;
//   const CFreal myGhost = 2.0*alpha*By0 - my;
//   const CFreal mzGhost = 2.0*alpha*Bz0 - mz;
//   const CFreal rhoGhost = _dataGhostState[MHDTerm::RHO];

//   _dataGhostState[MHDTerm::VX] = mxGhost/rhoGhost;
//   _dataGhostState[MHDTerm::VY] = myGhost/rhoGhost;
//   _dataGhostState[MHDTerm::VZ] = mzGhost/rhoGhost;
//   _dataGhostState[MHDTerm::V] = sqrt(_dataGhostState[MHDTerm::VX]*_dataGhostState[MHDTerm::VX] +
//    _dataGhostState[MHDTerm::VY]*_dataGhostState[MHDTerm::VY] +
//    _dataGhostState[MHDTerm::VZ]*_dataGhostState[MHDTerm::VZ]);

  /////////////////////////////////////////////////////////////////////////////////

  const CFreal bn = _dataInnerState[MHDTerm::BX]*nx +
    _dataInnerState[MHDTerm::BY]*ny +
    _dataInnerState[MHDTerm::BZ]*nz;

  RealVector BDipoleInner(PhysicalModelStack::getActive()->getDim()), 
    BDipoleGhost(PhysicalModelStack::getActive()->getDim()),
    ghostStateCoord(PhysicalModelStack::getActive()->getDim());

  ghostStateCoord = ghostState->getCoordinates();
  const CFreal mX = _varSet->getMX();
  const CFreal mY = _varSet->getMY();
  const CFreal mZ = _varSet->getMZ();

  CFreal vXGhost = 0.;
  CFreal vYGhost = 0.;
  CFreal vZGhost = 0.;

  const CFreal pi = 3.141592653;

  // e.g. period of rotation of earth (1 day (in seconds) non-dimensionalized with 15.9325 s (found by 1 earth radius (in km) / 400km/s))
  // _T = 86400 / 15.9325;

  if ((mX == 0.0) && (mY == 0.0)) {
    vXGhost = 2.0*pi*ghostStateCoord[1]/_T;
    vYGhost = -2.0*pi*ghostStateCoord[0]/_T;
    vZGhost = 0.0;
  } 
  else if ((mX != 0.0) && (mY == 0.0)) {
         const CFreal theta = atan(mZ/mX);
         vXGhost = 2.0*pi*ghostStateCoord[1]*sin(theta)/_T;
         vYGhost = 2.0*pi*(ghostStateCoord[2]*cos(theta)-ghostStateCoord[0]*sin(theta))/_T;
         vZGhost = 2.0*pi*ghostStateCoord[1]*cos(theta)/_T;
       } 
       else if ((mX != 0.0) && (mY != 0.0)) {
              const CFreal theta = atan(mZ/sqrt(mX*mX+mY*mY));
              const CFreal phi = atan(mX/mY);
              vXGhost = 2.0*pi*(ghostStateCoord[1]*sin(theta)-ghostStateCoord[2]*cos(theta)*cos(phi))/_T;
              vYGhost = 2.0*pi*(ghostStateCoord[2]*cos(theta)*sin(phi)-ghostStateCoord[0]*sin(theta))/_T;
              vZGhost = 2.0*pi*(ghostStateCoord[0]*cos(theta)*cos(phi)-ghostStateCoord[1]*sin(phi)*cos(theta))/_T;
       }
  
  const CFreal vGhost = sqrt(vXGhost*vXGhost + vYGhost*vYGhost + vZGhost*vZGhost);    

  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);

  // all boundary conditions are according to Gabor Toth's suggestions
  _dataGhostState[MHDTerm::RHO] = _rhoFixed;
  _dataGhostState[MHDTerm::VX] = vXGhost;
  _dataGhostState[MHDTerm::VY] = vYGhost;
  _dataGhostState[MHDTerm::VZ] = vZGhost;
  _dataGhostState[MHDTerm::BX] = _dataInnerState[MHDTerm::BX] - 2.0*bn*nx;
  _dataGhostState[MHDTerm::BY] = _dataInnerState[MHDTerm::BY] - 2.0*bn*ny;
  _dataGhostState[MHDTerm::BZ] = _dataInnerState[MHDTerm::BZ] - 2.0*bn*nz; 
  _dataGhostState[MHDTerm::V] = vGhost;
  _dataGhostState[MHDTerm::P] = _dataInnerState[MHDTerm::P];
  _dataGhostState[MHDTerm::A] = sqrt(_varSet->getModel()->getGamma()*_dataGhostState[MHDTerm::P]/_dataGhostState[MHDTerm::RHO]);
  _dataGhostState[MHDTerm::B] = _dataInnerState[MHDTerm::B];
  _dataGhostState[MHDProjectionTerm::PHI] = _dataInnerState[MHDProjectionTerm::PHI];
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
