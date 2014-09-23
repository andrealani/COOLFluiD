#include "FiniteVolumeMaxwell/FiniteVolumeMaxwell.hh"
#include "TruncatedDomain2DProjection.hh"
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

MethodCommandProvider<TruncatedDomain2DProjection, CellCenterFVMData, FiniteVolumeMaxwellModule> truncatedDomain2DProjectionFVMCCProvider("TruncatedDomain2DProjectionFVMCC");

//////////////////////////////////////////////////////////////////////////////

TruncatedDomain2DProjection::TruncatedDomain2DProjection(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
}

//////////////////////////////////////////////////////////////////////////////

TruncatedDomain2DProjection::~TruncatedDomain2DProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

void TruncatedDomain2DProjection::setup()
{
  FVMCC_BC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<Maxwell2DProjectionAdimVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void TruncatedDomain2DProjection::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
//   const CFuint faceID = face->getID();
//   const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
   
//    DataHandle<CFreal> normals = socket_normals.getDataHandle();
//    
//    CFreal nx = normals[startID];
//    CFreal ny = normals[startID + 1];
//    CFreal nz = 0;

  cf_assert(_varSet.isNotNull());
  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  _dataGhostState[ConvMaxwellTerm::BX] = _dataInnerState[ConvMaxwellTerm::BX];
  _dataGhostState[ConvMaxwellTerm::BY] = _dataInnerState[ConvMaxwellTerm::BY];
  _dataGhostState[ConvMaxwellTerm::BZ] = _dataInnerState[ConvMaxwellTerm::BZ];
  _dataGhostState[ConvMaxwellTerm::EX] = _dataInnerState[ConvMaxwellTerm::EX];
  _dataGhostState[ConvMaxwellTerm::EY] = _dataInnerState[ConvMaxwellTerm::EY];
  _dataGhostState[ConvMaxwellTerm::EZ] = _dataInnerState[ConvMaxwellTerm::EZ];
  _dataGhostState[MaxwellProjectionAdimTerm::PSI] = -_dataInnerState[MaxwellProjectionAdimTerm::PSI];	//AS done in Munz2
  _dataGhostState[MaxwellProjectionAdimTerm::PHI] = -_dataInnerState[MaxwellProjectionAdimTerm::PHI];	//As done in Munz2
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
