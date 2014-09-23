#include "FiniteVolumeMaxwell/FiniteVolumeMaxwell.hh"
#include "PerfectConductingWall3DProjectionDim.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Maxwell;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PerfectConductingWall3DProjectionDim, CellCenterFVMData, FiniteVolumeMaxwellModule> perfectConductingWall3DProjectionDimFVMCCProvider("PerfectConductingWall3DProjectionDimFVMCC");

//////////////////////////////////////////////////////////////////////////////

PerfectConductingWall3DProjectionDim::PerfectConductingWall3DProjectionDim(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
}

//////////////////////////////////////////////////////////////////////////////

PerfectConductingWall3DProjectionDim::~PerfectConductingWall3DProjectionDim()
{
}

//////////////////////////////////////////////////////////////////////////////

void PerfectConductingWall3DProjectionDim::setup()
{
  FVMCC_BC::setup();
  
//   std::cout << "PerfectConductingWall3DProjectionDim::setup() \n";  
  _varSet = getMethodData().getUpdateVar().d_castTo<Maxwell3DProjectionVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void PerfectConductingWall3DProjectionDim::setGhostState(GeometricEntity *const face)
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

  cf_assert(_varSet.isNotNull());
  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);
  
  const CFreal bn = _dataInnerState[ConvMaxwellTerm::BX]*nx +
    _dataInnerState[ConvMaxwellTerm::BY]*ny + _dataInnerState[ConvMaxwellTerm::BZ]*nz;
  const CFreal en = _dataInnerState[ConvMaxwellTerm::EX]*nx +
    _dataInnerState[ConvMaxwellTerm::EY]*ny + _dataInnerState[ConvMaxwellTerm::EZ]*nz ;  
//   const CFreal chi = _varSet->getModel()->getDivECleaningConst();
//   const CFreal c_e = _varSet->getModel()->getLightSpeed();//speed of light
  


///Munz1
  _dataGhostState[ConvMaxwellTerm::BX] = _dataInnerState[ConvMaxwellTerm::BX] - 2*bn*nx;
  _dataGhostState[ConvMaxwellTerm::BY] = _dataInnerState[ConvMaxwellTerm::BY] - 2*bn*ny;
  _dataGhostState[ConvMaxwellTerm::BZ] = _dataInnerState[ConvMaxwellTerm::BZ] - 2*bn*nz;
  _dataGhostState[ConvMaxwellTerm::EX] = -_dataInnerState[ConvMaxwellTerm::EX] + 2*en*nx;
  _dataGhostState[ConvMaxwellTerm::EY] = -_dataInnerState[ConvMaxwellTerm::EY] + 2*en*ny;
  _dataGhostState[ConvMaxwellTerm::EZ] = -_dataInnerState[ConvMaxwellTerm::EZ] + 2*en*nz;
  _dataGhostState[MaxwellProjectionTerm::PSI] = _dataInnerState[MaxwellProjectionTerm::PSI];
  _dataGhostState[MaxwellProjectionTerm::PHI] = -_dataInnerState[MaxwellProjectionTerm::PHI];
     
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
