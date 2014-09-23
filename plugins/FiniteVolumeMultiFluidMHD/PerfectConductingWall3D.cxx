#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/PerfectConductingWall3D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MultiFluidMHD/DiffMFMHD3DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "Framework/MeshData.hh"
  
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace COOLFluiD::Physics::Maxwell;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PerfectConductingWall3D, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> perfectConductingWall3DFVMCCProvider("PerfectConductingWall3DFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void PerfectConductingWall3D::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

PerfectConductingWall3D::PerfectConductingWall3D
(const std::string& name) :
  FVMCC_BC(name),
  _updateVarSet(CFNULL)
{
}
      
//////////////////////////////////////////////////////////////////////////////

PerfectConductingWall3D::~PerfectConductingWall3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void PerfectConductingWall3D::setup()
{
  FVMCC_BC::setup();
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >();
}
      
//////////////////////////////////////////////////////////////////////////////

void PerfectConductingWall3D::setGhostState(GeometricEntity *const face)
{
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);
  
  State *const innerState = face->getState(0); 
  State *const ghostState = face->getState(1);
  
  
  ///Maxwell Equations Perfectly Conducting Wall Condition
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

  cf_assert(_updateVarSet.isNotNull());
  
  const CFreal bn = (*innerState)[0]*nx + (*innerState)[1]*ny + (*innerState)[2]*nz;
  const CFreal en = (*innerState)[3]*nx + (*innerState)[4]*ny + (*innerState)[5]*nz;  
//  const CFreal chi = _varSet->getModel()->getDivECleaningConst();

  (*ghostState)[0] = (*innerState)[0] - 2*bn*nx;	//Bx
  (*ghostState)[1] = (*innerState)[1] - 2*bn*ny;	//By
  (*ghostState)[2] = (*innerState)[2] - 2*bn*nz;	//Bz
  (*ghostState)[3] = -(*innerState)[3] + 2*en*nx;	//Ex
  (*ghostState)[4] = -(*innerState)[4] + 2*en*ny;	//Ey
  (*ghostState)[5] = -(*innerState)[5] + 2*en*nz;	//Ez
  (*ghostState)[6] = (*innerState)[6];			//Psi
  (*ghostState)[7] = -(*innerState)[7];			//Phi
  
///MultiFluidMHD mirror Condition in 3D
  const CFuint endEM = 8;
 
  //set the densities
  for (CFuint i = 0 ; i < nbSpecies; i++){
    (*ghostState)[endEM + i] = (*innerState)[endEM + i];
  }
 
  //set the Velocities
  CFuint dim = 3;
  for (CFuint i = 0 ; i < nbSpecies; i++){
    CFreal un_i = (*innerState)[endEM + nbSpecies + dim*i]*nx + (*innerState)[endEM + nbSpecies + dim*i + 1]*ny + (*innerState)[endEM + nbSpecies + dim*i + 2]*nz; 
    (*ghostState)[endEM + nbSpecies + dim*i] = (*innerState)[endEM + nbSpecies + dim*i] - 2*un_i*nx;
    (*ghostState)[endEM + nbSpecies + dim*i + 1] = (*innerState)[endEM + nbSpecies + dim*i + 1] - 2*un_i*ny;
    (*ghostState)[endEM + nbSpecies + dim*i + 2] = (*innerState)[endEM + nbSpecies + dim*i + 2] - 2*un_i*nz;
  } 
 
  //set the Temperatures
  for (CFuint i = 0 ; i < nbSpecies; i++){
    (*ghostState)[endEM + nbSpecies + dim*nbSpecies + i] = (*innerState)[endEM + nbSpecies + dim*nbSpecies + i];
    cf_assert((*innerState)[endEM + nbSpecies + dim*nbSpecies + i] > 0.); 
  } 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
