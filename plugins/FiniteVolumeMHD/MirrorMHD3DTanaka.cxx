#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MirrorMHD3DTanaka.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD3DVarSet.hh"
#include "Framework/MeshData.hh"
  
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

MethodCommandProvider<MirrorMHD3DTanaka, CellCenterFVMData, FiniteVolumeMHDModule> 
mirrorMHD3DTanakaFVMCCProvider("MirrorMHD3DTanakaFVMCC");

//////////////////////////////////////////////////////////////////////
   
void MirrorMHD3DTanaka::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("rhoFixed","Density value that is to be fixed when the normal velocity is outward from the body.");
}

//////////////////////////////////////////////////////////////////////////////

MirrorMHD3DTanaka::MirrorMHD3DTanaka(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
  addConfigOptionsTo(this);

  _rhoFixed = 1.0;
  setParameter("rhoFixed",&_rhoFixed);
}

//////////////////////////////////////////////////////////////////////////////

MirrorMHD3DTanaka::~MirrorMHD3DTanaka()
{
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DTanaka::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void MirrorMHD3DTanaka::setup()
{
  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void MirrorMHD3DTanaka::setGhostState(GeometricEntity *const face)
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
    _dataInnerState[MHDTerm::BY]*ny + _dataInnerState[MHDTerm::BZ]*nz;

  
  // all boundary conditions are according to Gabor Toth's suggestions
  _dataGhostState[MHDTerm::RHO] = _rhoFixed;
  _dataGhostState[MHDTerm::VX] = -_dataInnerState[MHDTerm::VX];
  _dataGhostState[MHDTerm::VY] = -_dataInnerState[MHDTerm::VY];
  _dataGhostState[MHDTerm::VZ] = -_dataInnerState[MHDTerm::VZ];
  _dataGhostState[MHDTerm::BX] = _dataInnerState[MHDTerm::BX] - 2.0*bn*nx;
  _dataGhostState[MHDTerm::BY] = _dataInnerState[MHDTerm::BY] - 2.0*bn*ny;
  _dataGhostState[MHDTerm::BZ] = _dataInnerState[MHDTerm::BZ] - 2.0*bn*nz; 
  _dataGhostState[MHDTerm::V] = _dataInnerState[MHDTerm::V];
  _dataGhostState[MHDTerm::P] = _dataInnerState[MHDTerm::P];
  _dataGhostState[MHDTerm::A] = sqrt(_varSet->getModel()->getGamma()*_dataGhostState[MHDTerm::P]/_dataGhostState[MHDTerm::RHO]);
  _dataGhostState[MHDTerm::B] = _dataInnerState[MHDTerm::B];
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
