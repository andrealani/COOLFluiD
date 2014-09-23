#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MirrorMHD3DTanakaPFix.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD3DVarSet.hh"
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

MethodCommandProvider<MirrorMHD3DTanakaPFix, CellCenterFVMData, FiniteVolumeMHDModule> 
mirrorMHD3DTanakaPFixFVMCCProvider("MirrorMHD3DTanakaPFixFVMCC");

//////////////////////////////////////////////////////////////////////
   
void MirrorMHD3DTanakaPFix::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >
     ("rhoFixed", "Density value that is to be fixed for introducing numerical diffusion to get non-zero r-momentum thus imitating the plasma flow at the ionosphere/magnetosphere interface.");
  options.addConfigOption< CFreal >
     ("pFixed", "Pressure value that is to be fixed according to Powell.");
}

//////////////////////////////////////////////////////////////////////

MirrorMHD3DTanakaPFix::MirrorMHD3DTanakaPFix(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
  addConfigOptionsTo(this);

  _rhoFixed = 1.0;
  setParameter("rhoFixed",&_rhoFixed);

  _pFixed = 1.0;
  setParameter("pFixed",&_pFixed);
}

//////////////////////////////////////////////////////////////////////

MirrorMHD3DTanakaPFix::~MirrorMHD3DTanakaPFix()
{
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DTanakaPFix::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DTanakaPFix::setup()
{
 FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DTanakaPFix::setGhostState(GeometricEntity *const face)
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

  const CFreal bn = _dataInnerState[MHDTerm::BX]*nx +
    _dataInnerState[MHDTerm::BY]*ny +
    _dataInnerState[MHDTerm::BZ]*nz;

  // all boundary conditions are according to Gabor Toth's (and for pressure fixation Powell's) suggestions
  
  _dataGhostState[MHDTerm::RHO] = _rhoFixed;
  _dataGhostState[MHDTerm::VX] = -_dataInnerState[MHDTerm::VX];
  _dataGhostState[MHDTerm::VY] = -_dataInnerState[MHDTerm::VY];
  _dataGhostState[MHDTerm::VZ] = -_dataInnerState[MHDTerm::VZ];
  _dataGhostState[MHDTerm::BX] = _dataInnerState[MHDTerm::BX] - 2.0*bn*nx;
  _dataGhostState[MHDTerm::BY] = _dataInnerState[MHDTerm::BY] - 2.0*bn*ny;
  _dataGhostState[MHDTerm::BZ] = _dataInnerState[MHDTerm::BZ] - 2.0*bn*nz; 
  _dataGhostState[MHDTerm::V] = _dataInnerState[MHDTerm::V];
  _dataGhostState[MHDTerm::P] = _pFixed;
  _dataGhostState[MHDTerm::A] = sqrt(_varSet->getModel()->getGamma()*_dataGhostState[MHDTerm::P]/_dataGhostState[MHDTerm::RHO]);
  _dataGhostState[MHDTerm::B] = _dataInnerState[MHDTerm::B];
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
