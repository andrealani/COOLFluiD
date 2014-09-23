#include "FiniteVolumeMaxwell/FiniteVolumeMaxwell.hh"
#include "UnsteadyPerfectConductingWall2DProjectionDim.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
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

MethodCommandProvider<UnsteadyPerfectConductingWall2DProjectionDim, CellCenterFVMData, FiniteVolumeMaxwellModule> unsteadyPerfectConductingWall2DProjectionDimFVMCCProvider("UnsteadyPerfectConductingWall2DProjectionDimFVMCC");

//////////////////////////////////////////////////////////////////////////////

UnsteadyPerfectConductingWall2DProjectionDim::UnsteadyPerfectConductingWall2DProjectionDim(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyPerfectConductingWall2DProjectionDim::~UnsteadyPerfectConductingWall2DProjectionDim()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyPerfectConductingWall2DProjectionDim::setup()
{
  FVMCC_BC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<Maxwell2DProjectionVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyPerfectConductingWall2DProjectionDim::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
   const CFuint faceID = face->getID();
   const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
   
   DataHandle<CFreal> normals = socket_normals.getDataHandle();
   
   CFreal nx = normals[startID];
   CFreal ny = normals[startID + 1];
   CFreal nz = 0;
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx *= invFaceLength;
  ny *= invFaceLength;
  nz *= invFaceLength;

  cf_assert(_varSet.isNotNull());
  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);
  
  const CFreal bn = _dataInnerState[ConvMaxwellTerm::BX]*nx +
    _dataInnerState[ConvMaxwellTerm::BY]*ny;
  const CFreal en = _dataInnerState[ConvMaxwellTerm::EX]*nx +
    _dataInnerState[ConvMaxwellTerm::EY]*ny;  
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
