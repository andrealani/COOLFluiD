#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MirrorMHD2DProjection.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD2DProjectionVarSet.hh"
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

MethodCommandProvider<MirrorMHD2DProjection, CellCenterFVMData, FiniteVolumeMHDModule> 
mirrorMHD2DProjectionFVMCCProvider("MirrorMHD2DProjectionFVMCC");

//////////////////////////////////////////////////////////////////////////////

MirrorMHD2DProjection::MirrorMHD2DProjection(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
}

//////////////////////////////////////////////////////////////////////////////

MirrorMHD2DProjection::~MirrorMHD2DProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

void MirrorMHD2DProjection::setup()
{
  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD2DProjectionVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void MirrorMHD2DProjection::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
  
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);
  
  const CFreal vn = _dataInnerState[MHDTerm::VX]*nx + 
    _dataInnerState[MHDTerm::VY]*ny;
  const CFreal bn = _dataInnerState[MHDTerm::BX]*nx +
    _dataInnerState[MHDTerm::BY]*ny;
  
  _dataGhostState[MHDTerm::RHO] = _dataInnerState[MHDTerm::RHO];
  _dataGhostState[MHDTerm::VX] = _dataInnerState[MHDTerm::VX] - 2.0*vn*nx;
  _dataGhostState[MHDTerm::VY] = _dataInnerState[MHDTerm::VY] - 2.0*vn*ny;
  _dataGhostState[MHDTerm::VZ] = _dataInnerState[MHDTerm::VZ];
  _dataGhostState[MHDTerm::BX] = _dataInnerState[MHDTerm::BX] - 2.0*bn*nx;
  _dataGhostState[MHDTerm::BY] = _dataInnerState[MHDTerm::BY] - 2.0*bn*ny;
  _dataGhostState[MHDTerm::BZ] = _dataInnerState[MHDTerm::BZ]; 
  _dataGhostState[MHDTerm::V] = _dataInnerState[MHDTerm::V];
  _dataGhostState[MHDTerm::A] = _dataInnerState[MHDTerm::A];
  _dataGhostState[MHDTerm::P] = _dataInnerState[MHDTerm::P];
  _dataGhostState[MHDTerm::B] = _dataInnerState[MHDTerm::B];
  _dataGhostState[MHDProjectionTerm::PHI] = _dataInnerState[MHDProjectionTerm::PHI];
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
