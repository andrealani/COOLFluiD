#include "FiniteVolumeMaxwell/FiniteVolumeMaxwell.hh"
#include "UnsteadyPerfectConductingWall2DProjection.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Maxwell/Maxwell2DProjectionAdimVarSet.hh"
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

MethodCommandProvider<UnsteadyPerfectConductingWall2DProjection, CellCenterFVMData, FiniteVolumeMaxwellModule> unsteadyPerfectConductingWall2DProjectionFVMCCProvider("UnsteadyPerfectConductingWall2DProjectionFVMCC");

//////////////////////////////////////////////////////////////////////////////

UnsteadyPerfectConductingWall2DProjection::UnsteadyPerfectConductingWall2DProjection(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyPerfectConductingWall2DProjection::~UnsteadyPerfectConductingWall2DProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyPerfectConductingWall2DProjection::setup()
{
  FVMCC_BC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<Maxwell2DProjectionAdimVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyPerfectConductingWall2DProjection::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
   const CFuint faceID = face->getID();
   const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
   
   DataHandle<CFreal> normals = socket_normals.getDataHandle();
   
   CFreal nx = normals[startID];
   CFreal ny = normals[startID + 1];
   CFreal nz = 0;
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;

  cf_assert(_varSet.isNotNull());
  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);
  
  const CFreal bn = _dataInnerState[ConvMaxwellTerm::BX]*nx +
    _dataInnerState[ConvMaxwellTerm::BY]*ny;
  const CFreal en = _dataInnerState[ConvMaxwellTerm::EX]*nx +
    _dataInnerState[ConvMaxwellTerm::EY]*ny;  
  const CFreal chi = _varSet->getModel()->getDivECleaningConst();


///Munz1
  _dataGhostState[ConvMaxwellTerm::BX] = _dataInnerState[ConvMaxwellTerm::BX] - 2*bn*nx;
  _dataGhostState[ConvMaxwellTerm::BY] = _dataInnerState[ConvMaxwellTerm::BY] - 2*bn*ny;
  _dataGhostState[ConvMaxwellTerm::BZ] = _dataInnerState[ConvMaxwellTerm::BZ] - 2*bn*nz;
  _dataGhostState[ConvMaxwellTerm::EX] = -_dataInnerState[ConvMaxwellTerm::EX] + 2*en*nx;
  _dataGhostState[ConvMaxwellTerm::EY] = -_dataInnerState[ConvMaxwellTerm::EY] + 2*en*ny;
  _dataGhostState[ConvMaxwellTerm::EZ] = -_dataInnerState[ConvMaxwellTerm::EZ] + 2*en*nz;
  _dataGhostState[MaxwellProjectionAdimTerm::PSI] = _dataInnerState[MaxwellProjectionAdimTerm::PSI];
  _dataGhostState[MaxwellProjectionAdimTerm::PHI] = -_dataInnerState[MaxwellProjectionAdimTerm::PHI];
    
///Mun2 transmission case
  
//   _dataGhostState[ConvMaxwellTerm::BX] = _dataInnerState[ConvMaxwellTerm::BX] - 2*bn*nx;
//   _dataGhostState[ConvMaxwellTerm::BY] = _dataInnerState[ConvMaxwellTerm::BY] - 2*bn*ny;
//   _dataGhostState[ConvMaxwellTerm::BZ] = _dataInnerState[ConvMaxwellTerm::BZ];
//   _dataGhostState[ConvMaxwellTerm::EX] = -_dataInnerState[ConvMaxwellTerm::EX];
//   _dataGhostState[ConvMaxwellTerm::EY] = -_dataInnerState[ConvMaxwellTerm::EY];
//   _dataGhostState[ConvMaxwellTerm::EZ] = -_dataInnerState[ConvMaxwellTerm::EZ];
//   _dataGhostState[MaxwellProjectionAdimTerm::PSI] = -_dataInnerState[MaxwellProjectionAdimTerm::PSI] ;
//   _dataGhostState[MaxwellProjectionAdimTerm::PHI] = _dataInnerState[MaxwellProjectionAdimTerm::PHI] - 2*chi*(nx*_dataInnerState[ConvMaxwellTerm::EX] + ny*_dataInnerState[ConvMaxwellTerm::EY]);
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
