#include "FiniteVolumeMaxwell/FiniteVolumeMaxwell.hh"
#include "SimmetryPlane2DProjectionDim.hh"
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

MethodCommandProvider<SimmetryPlane2DProjectionDim, CellCenterFVMData, FiniteVolumeMaxwellModule> simmetryPlane2DProjectionDimFVMCCProvider("SimmetryPlane2DProjectionDimFVMCC");

///This options could be added in the future if it works
// //////////////////////////////////////////////////////////////////////////////
// 
// void SimmetryPlane2DProjectionDim::defineConfigOptions(Config::OptionList& options)
// {
//   options.addConfigOption< CFreal >
//     ("Ez0","Imposed Ez0");
// 
//   options.addConfigOption< std::vector<CFreal> >
//     ("electricalConductivity","electrical Conductivity used in Ohm's law");
// 
// }

//////////////////////////////////////////////////////////////////////////////

SimmetryPlane2DProjectionDim::SimmetryPlane2DProjectionDim(const std::string& name) :
  FVMCC_BC(name)/*,
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()*/
{
//   addConfigOptionsTo(this);
//   
//   _electricalConductivity = 66938.88;
//   setParameter("electricalConductivity",&_electricalConductivity);
//   
//   _Ez0 = -0.1;
//   setParameter("Ez0",&_Ez0);
// 
//  
}

//////////////////////////////////////////////////////////////////////////////

SimmetryPlane2DProjectionDim::~SimmetryPlane2DProjectionDim()
{
}

//////////////////////////////////////////////////////////////////////////////

void SimmetryPlane2DProjectionDim::setup()
{
  FVMCC_BC::setup();
}

//////////////////////////////////////////////////////////////////////////////

void SimmetryPlane2DProjectionDim::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
//   const CFreal muZero = 1.2566370614359e-6;
//   const CFreal omega = 66938.88;
//   const CFreal Ez0 = -0.1;  
 
  
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = 0;
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx *= invFaceLength;
  ny *= invFaceLength;
  nz *= invFaceLength;

  // here refactor to use user-defined indices, one for the dB[idxB]/dx[idxX] 
  const CFreal dy = ghostState->getCoordinates()[YY] - innerState->getCoordinates()[YY];
  
///Munz1
  (*ghostState)[ConvMaxwellTerm::BX] = -(*innerState)[ConvMaxwellTerm::BX];
  (*ghostState)[ConvMaxwellTerm::BY] = (*innerState)[ConvMaxwellTerm::BY] ;
  (*ghostState)[ConvMaxwellTerm::BZ] = -(*innerState)[ConvMaxwellTerm::BZ];
  (*ghostState)[ConvMaxwellTerm::EX] = (*innerState)[ConvMaxwellTerm::EX];
  (*ghostState)[ConvMaxwellTerm::EY] = -(*innerState)[ConvMaxwellTerm::EY];
  (*ghostState)[ConvMaxwellTerm::EZ] = (*innerState)[ConvMaxwellTerm::EZ];
  (*ghostState)[MaxwellProjectionAdimTerm::PSI] = (*innerState)[MaxwellProjectionAdimTerm::PSI];
  (*ghostState)[MaxwellProjectionAdimTerm::PHI] = (*innerState)[MaxwellProjectionAdimTerm::PHI];
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
