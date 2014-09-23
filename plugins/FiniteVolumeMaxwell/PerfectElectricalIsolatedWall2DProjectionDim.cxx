#include "FiniteVolumeMaxwell/FiniteVolumeMaxwell.hh"
#include "PerfectElectricalIsolatedWall2DProjectionDim.hh"
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

MethodCommandProvider<PerfectElectricalIsolatedWall2DProjectionDim, CellCenterFVMData, FiniteVolumeMaxwellModule> perfectElectricalIsolatedWall2DProjectionDimFVMCCProvider("PerfectElectricalIsolatedWall2DProjectionDimFVMCC");

//////////////////////////////////////////////////////////////////////////////

PerfectElectricalIsolatedWall2DProjectionDim::PerfectElectricalIsolatedWall2DProjectionDim(const std::string& name) :
  FVMCC_BC(name)/*,
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()*/
{
}

//////////////////////////////////////////////////////////////////////////////

PerfectElectricalIsolatedWall2DProjectionDim::~PerfectElectricalIsolatedWall2DProjectionDim()
{
}

//////////////////////////////////////////////////////////////////////////////

void PerfectElectricalIsolatedWall2DProjectionDim::setup()
{
  FVMCC_BC::setup();
  
//   _varSet = getMethodData().getUpdateVar().d_castTo<Maxwell2DProjectionVarSet>();
//   _varSet->getModel()->resizePhysicalData(_dataInnerState);
//   _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void PerfectElectricalIsolatedWall2DProjectionDim::setGhostState(GeometricEntity *const face)
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

//   cf_assert(_varSet.isNotNull());
  // set the physical data starting from the inner state
//   _varSet->computePhysicalData(*innerState, _dataInnerState);
  
  const CFreal bn = (*innerState)[ConvMaxwellTerm::BX]*nx +
    (*innerState)[ConvMaxwellTerm::BY]*ny;
/*  const CFreal en = _dataInnerState[ConvMaxwellTerm::EX]*nx +
    _dataInnerState[ConvMaxwellTerm::EY]*ny; */ 
//   const CFreal chi = _varSet->getModel()->getDivECleaningConst();
//   const CFreal c_e = _varSet->getModel()->getLightSpeed();//speed of light
  


///Munz1
  (*ghostState)[ConvMaxwellTerm::BX] = (*innerState)[ConvMaxwellTerm::BX] /*+ 2*bn*nx*/;
  (*ghostState)[ConvMaxwellTerm::BY] = -(*innerState)[ConvMaxwellTerm::BY] /*+ 2*bn*ny*/;
  (*ghostState)[ConvMaxwellTerm::BZ] = -(*innerState)[ConvMaxwellTerm::BZ] + 2*bn*nz;
  (*ghostState)[ConvMaxwellTerm::EX] = -(*innerState)[ConvMaxwellTerm::EX]/* + 2*en*nx*/;
  (*ghostState)[ConvMaxwellTerm::EY] = -(*innerState)[ConvMaxwellTerm::EY]/* + 2*en*ny*/;
  (*ghostState)[ConvMaxwellTerm::EZ] = -(*innerState)[ConvMaxwellTerm::EZ]/* + 2*en*nz*/;
  (*ghostState)[MaxwellProjectionTerm::PSI] = (*innerState)[MaxwellProjectionTerm::PSI];
  (*ghostState)[MaxwellProjectionTerm::PHI] = -(*innerState)[MaxwellProjectionTerm::PHI];
     
//   _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
