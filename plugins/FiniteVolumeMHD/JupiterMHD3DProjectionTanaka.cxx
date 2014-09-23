#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "JupiterMHD3DProjectionTanaka.hh"
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

MethodCommandProvider<JupiterMHD3DProjectionTanaka, CellCenterFVMData, FiniteVolumeMHDModule> jupiterMHD3DProjectionTanakaFVMCCProvider("JupiterMHD3DProjectionTanakaFVMCC");

//////////////////////////////////////////////////////////////////////
   
void JupiterMHD3DProjectionTanaka::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("rhoFixed","Density value that is to be fixed when the normal velocity is outward from the body.");
}

//////////////////////////////////////////////////////////////////////

JupiterMHD3DProjectionTanaka::JupiterMHD3DProjectionTanaka(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
  addConfigOptionsTo(this);

  _rhoFixed = 1.0;
  setParameter("rhoFixed",&_rhoFixed);
}

//////////////////////////////////////////////////////////////////////

JupiterMHD3DProjectionTanaka::~JupiterMHD3DProjectionTanaka()
{
}

//////////////////////////////////////////////////////////////////////

void JupiterMHD3DProjectionTanaka::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////

void JupiterMHD3DProjectionTanaka::setup()
{
  FVMCC_BC::setup();
  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////

void JupiterMHD3DProjectionTanaka::setGhostState(GeometricEntity *const face)
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

  const CFreal vn = _dataInnerState[MHDTerm::VX]*nx +
    _dataInnerState[MHDTerm::VY]*ny +
    _dataInnerState[MHDTerm::VZ]*nz;

  _dataGhostState[MHDTerm::RHO] = _rhoFixed;
  _dataGhostState[MHDTerm::VX] = _dataInnerState[MHDTerm::VX] - 2.0*vn*nx;
  _dataGhostState[MHDTerm::VY] = _dataInnerState[MHDTerm::VY] - 2.0*vn*ny;
  _dataGhostState[MHDTerm::VZ] = _dataInnerState[MHDTerm::VZ] - 2.0*vn*nz;
  _dataGhostState[MHDTerm::BX] = _dataInnerState[MHDTerm::BX];
  _dataGhostState[MHDTerm::BY] = _dataInnerState[MHDTerm::BY];
  _dataGhostState[MHDTerm::BZ] = _dataInnerState[MHDTerm::BZ]; 
  _dataGhostState[MHDTerm::V] = _dataInnerState[MHDTerm::V];
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
