#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/SubOutletPPCWRhoiViTi3D.hh"
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

MethodCommandProvider<SubOutletPPCWRhoiViTi3D, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> subOutletPPCWRhoiViTi3DFVMCCProvider("SubOutletPPCWRhoiViTi3DFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void SubOutletPPCWRhoiViTi3D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<CFreal> >("Pi","static pressure of the species");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletPPCWRhoiViTi3D::SubOutletPPCWRhoiViTi3D
(const std::string& name) :
  FVMCC_BC(name),
  _updateVarSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
   _Pi = std::vector<CFreal>();
   setParameter("Pi",&_Pi);
}

//////////////////////////////////////////////////////////////////////////////

SubOutletPPCWRhoiViTi3D::~SubOutletPPCWRhoiViTi3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletPPCWRhoiViTi3D::setup()
{
  FVMCC_BC::setup();
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >();
  _updateVarSet->getModel()->resizePhysicalData(_dataInnerState);
  _updateVarSet->getModel()->resizePhysicalData(_dataGhostState);
  
//   const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);  
//   for (CFuint i = 0; i < nbSpecies; ++i) {
//     cf_assert(_wallTemp[i] >= 0.0);
//   }

  // adimensionalize the temperature
//   _wallTemp /= _updateVarSet->getModel()->getTempRef();
}
      
//////////////////////////////////////////////////////////////////////////////

void SubOutletPPCWRhoiViTi3D::setGhostState(GeometricEntity *const face)
{
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);
  
  State *const innerState = face->getState(0); 
  State *const ghostState = face->getState(1);
    
  // set the physical data starting from the inner state
  _updateVarSet->computePhysicalData(*innerState, _dataInnerState);
  
  
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

  (_dataGhostState)[0] = (_dataInnerState)[0] - 2*bn*nx;	//Bx
  (_dataGhostState)[1] = (_dataInnerState)[1] - 2*bn*ny;	//By
  (_dataGhostState)[2] = (_dataInnerState)[2] - 2*bn*nz;	//Bz
  (_dataGhostState)[3] = -(_dataInnerState)[3] + 2*en*nx;	//Ex
  (_dataGhostState)[4] = -(_dataInnerState)[4] + 2*en*ny;	//Ey
  (_dataGhostState)[5] = -(_dataInnerState)[5] + 2*en*nz;	//Ez
  (_dataGhostState)[6] = (_dataInnerState)[6];			//Psi
  (_dataGhostState)[7] = -(_dataInnerState)[7];			//Phi
  
///MultiFluidMHD Subsonic Outlet Condition imposing pressure in 3D
 // here a fix is needed in order to have always ghostT > 0
 // if ghostT < 0  then the inner value is set
  const CFuint endEM = 8;
  const CFuint firstSpecies = _updateVarSet->getModel()->getFirstScalarVar(0);  
  const CFuint firstVelocity = _updateVarSet->getModel()->getFirstScalarVar(1);   
  const CFuint firstTemperature = _updateVarSet->getModel()->getFirstScalarVar(2); 
  const CFreal gamma = _updateVarSet->getModel()->getGamma(); 
  const CFreal gammaDivGammaMinus1 = gamma/(gamma - 1.0);
  CFuint dim = 3;
 
  (_dataGhostState)[endEM] = (_dataInnerState)[endEM]; 							//RHO
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    const CFreal Pi = _Pi[ie];
    (_dataGhostState)[firstSpecies + ie] = (_dataInnerState)[firstSpecies + ie];			//yi
    (_dataGhostState)[firstVelocity + dim*ie] = (_dataInnerState)[firstVelocity + dim*ie];			//Vxi
    (_dataGhostState)[firstVelocity + dim*ie + 1] = (_dataInnerState)[firstVelocity + dim*ie + 1];		//Vyi
    (_dataGhostState)[firstVelocity + dim*ie + 2] = (_dataInnerState)[firstVelocity + dim*ie + 2];		//Vzi
    (_dataGhostState)[firstTemperature + 4*ie] = (_dataInnerState)[firstTemperature + 4*ie];		//Ti
    (_dataGhostState)[firstTemperature + 4*ie + 1] = 2*Pi - 
						    (_dataInnerState)[firstTemperature + 4*ie + 1];	//Pi
						    
    const CFreal Vi2 = (_dataGhostState)[firstVelocity + dim*ie]*(_dataGhostState)[firstVelocity + dim*ie] +
			(_dataGhostState)[firstVelocity + dim*ie + 1]*(_dataGhostState)[firstVelocity + dim*ie + 1] + 
            (_dataGhostState)[firstVelocity + dim*ie + 2]*(_dataGhostState)[firstVelocity + dim*ie + 2];
    const CFreal rhoi =(_dataGhostState)[endEM]*(_dataGhostState)[firstSpecies + ie];
    
    (_dataGhostState)[firstTemperature + 4*ie + 2] = sqrt(gamma*(_dataGhostState)[firstTemperature + 4*ie + 1]/rhoi);		//ai
    (_dataGhostState)[firstTemperature + 4*ie + 3] = (0.5*rhoi*Vi2 + 
						    gammaDivGammaMinus1*(_dataGhostState)[firstTemperature + 4*ie + 1])/rhoi;	//Hi
  }    
  
  _updateVarSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
