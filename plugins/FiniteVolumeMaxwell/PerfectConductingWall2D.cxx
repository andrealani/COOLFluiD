#include "FiniteVolumeMaxwell/FiniteVolumeMaxwell.hh"
#include "PerfectConductingWall2D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Maxwell/Maxwell2DVarSet.hh"
#ifdef CF_HAVE_CUDA
#include "Maxwell/Maxwell2DProjectionVarSetT.hh"
#endif
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

MethodCommandProvider<PerfectConductingWall2D, CellCenterFVMData, FiniteVolumeMaxwellModule> perfectConductingWall2DFVMCCProvider("PerfectConductingWall2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

PerfectConductingWall2D::PerfectConductingWall2D(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
}

//////////////////////////////////////////////////////////////////////////////

PerfectConductingWall2D::~PerfectConductingWall2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void PerfectConductingWall2D::setup()
{
  FVMCC_BC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<Maxwell2DVarSet>();  //Maxwell2DVarSet
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void PerfectConductingWall2D::setGhostState(GeometricEntity *const face)
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
//    std::cout<< "nx = "<< nx <<"\t ny = "<< ny <<"\n";

  cf_assert(_varSet.isNotNull());
  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);
  
  const CFreal bn = _dataInnerState[ConvMaxwellTerm::BX]*nx +
    _dataInnerState[ConvMaxwellTerm::BY]*ny;
  const CFreal en = _dataInnerState[ConvMaxwellTerm::EX]*nx +
    _dataInnerState[ConvMaxwellTerm::EY]*ny;    

/// Munz1
  _dataGhostState[ConvMaxwellTerm::BX] = _dataInnerState[ConvMaxwellTerm::BX] - 2*bn*nx;
  _dataGhostState[ConvMaxwellTerm::BY] = _dataInnerState[ConvMaxwellTerm::BY] - 2*bn*ny;
  _dataGhostState[ConvMaxwellTerm::BZ] = _dataInnerState[ConvMaxwellTerm::BZ] - 2*bn*nz;
//    std::cout<< "BzR = "<< _dataGhostState[ConvMaxwellTerm::BZ] <<"\t BzL = "<< _dataInnerState[ConvMaxwellTerm::BZ] <<"\n";
  _dataGhostState[ConvMaxwellTerm::EX] = -_dataInnerState[ConvMaxwellTerm::EX] + 2*en*nx;
//    std::cout<< "ExR = "<< _dataGhostState[ConvMaxwellTerm::EX] <<"\t ExL = "<< _dataInnerState[ConvMaxwellTerm::EX] <<"\n";
  _dataGhostState[ConvMaxwellTerm::EY] = -_dataInnerState[ConvMaxwellTerm::EY] + 2*en*ny;
//    std::cout<< "EyR = "<< _dataGhostState[ConvMaxwellTerm::EY] <<"\t EyL = "<< _dataInnerState[ConvMaxwellTerm::EY] <<"\n";
  _dataGhostState[ConvMaxwellTerm::EZ] = -_dataInnerState[ConvMaxwellTerm::EZ] + 2*en*nz;
    
 
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
