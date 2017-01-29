#include "MathTools/MathChecks.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/BCPeriodicMFMHD.hh"
#include "Common/MPI/MPIStructDef.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Framework/MeshData.hh"

#include <iostream>
#include <map>
#include <algorithm>

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

MethodCommandProvider<BCPeriodicMFMHD, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> 
BCPeriodicMFMHDFVMCCProvider("BCPeriodicMFMHDFVMCC");

//////////////////////////////////////////////////////////////////////////////

void BCPeriodicMFMHD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<CFreal> >("Pi","static pressure of the species");
}

//////////////////////////////////////////////////////////////////////////////

BCPeriodicMFMHD::BCPeriodicMFMHD(const std::string& name):
  BCPeriodic(name),
  _updateVarSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
  addConfigOptionsTo(this);
  
  _Pi = std::vector<CFreal>();
  setParameter("Pi",&_Pi);  
}

//////////////////////////////////////////////////////////////////////////////

BCPeriodicMFMHD::~BCPeriodicMFMHD()
{
}

//////////////////////////////////////////////////////////////////////////////

void BCPeriodicMFMHD::setup()
{
  CFLog(VERBOSE, "BCPeriodicMFMHD::setup() => start\n");
  
  /* Setup parent class */
  BCPeriodic::setup();
    
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  _updateVarSet->getModel()->resizePhysicalData(_dataInnerState);
  _updateVarSet->getModel()->resizePhysicalData(_dataGhostState);
  
  CFLog(VERBOSE, "BCPeriodicMFMHD::setup() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void BCPeriodicMFMHD::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  State *const periodicState = computePeriodicState(face);
  
  // set the physical data starting from the inner state
  _updateVarSet->computePhysicalData(*innerState, _dataInnerState);
  
  //Periodic Boundary conditions for Maxwell equations
  const CFuint endEM = 8;
  for(CFuint h=0; h < endEM; h++){
    (_dataInnerState)[h] = (*periodicState)[h];
  }
//   ///Maxwell Equations Perfectly Conducting Wall Condition
//   const CFuint faceID = face->getID();
//   const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
//    
//   DataHandle<CFreal> normals = socket_normals.getDataHandle();
//   CFreal nx = normals[startID];
//   CFreal ny = normals[startID + 1];
//   CFreal nz = 0;
//   const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
//   nx *= invFaceLength;
//   ny *= invFaceLength;
// 
//   cf_assert(_updateVarSet.isNotNull());
//   
//   const CFreal bn = (_dataInnerState)[0]*nx + (_dataInnerState)[1]*ny;
//   const CFreal en = (_dataInnerState)[3]*nx + (_dataInnerState)[4]*ny;  
// //  const CFreal chi = _varSet->getModel()->getDivECleaningConst();
// 
//   (_dataGhostState)[0] = -(_dataInnerState)[0] + 2*bn*nx;	//Bx
//   (_dataGhostState)[1] = -(_dataInnerState)[1] /*+ 2*bn*ny*/;	//By
//   (_dataGhostState)[2] = -(_dataInnerState)[2] + 2*bn*nz;	//Bz
//   (_dataGhostState)[3] = -(_dataInnerState)[3] /*- 2*en*nx*/;	//Ex
//   (_dataGhostState)[4] = -(_dataInnerState)[4] /*- 2*en*ny*/;	//Ey
//   (_dataGhostState)[5] = -(_dataInnerState)[5] /*- 2*en*nz*/;	//Ez
//   (_dataGhostState)[6] = (_dataInnerState)[6];			//Psi
//   (_dataGhostState)[7] = -(_dataInnerState)[7];			//Phi
  
  //MultiFluidMHD Subsonic Outlet Condition imposing pressure in 2D
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);

  const CFuint firstSpecies = _updateVarSet->getModel()->getFirstScalarVar(0);  
  const CFuint firstVelocity = _updateVarSet->getModel()->getFirstScalarVar(1);   
  const CFuint firstTemperature = _updateVarSet->getModel()->getFirstScalarVar(2); 
  const CFreal gamma = _updateVarSet->getModel()->getGamma(); 
  const CFreal gammaDivGammaMinus1 = gamma/(gamma - 1.0);
 
  (_dataGhostState)[endEM] = (_dataInnerState)[endEM]; 							//RHO
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    const CFreal Pi = _Pi[ie];
//     cout << "SubOutletPEIWRhoiViTi::setGhostState => Pi = " << Pi << endl;
    (_dataGhostState)[firstSpecies + ie] = (_dataInnerState)[firstSpecies + ie];			//yi
    (_dataGhostState)[firstVelocity + 2*ie] = (_dataInnerState)[firstVelocity + 2*ie];			//Vxi
    (_dataGhostState)[firstVelocity + 2*ie + 1] = (_dataInnerState)[firstVelocity + 2*ie + 1];		//Vyi
    (_dataGhostState)[firstTemperature + 4*ie] = (_dataInnerState)[firstTemperature + 4*ie];		//Ti
    (_dataGhostState)[firstTemperature + 4*ie + 1] = 2*Pi - 
						    (_dataInnerState)[firstTemperature + 4*ie + 1];	//Pi
						    
    const CFreal Vi2 = (_dataGhostState)[firstVelocity + 2*ie]*(_dataGhostState)[firstVelocity + 2*ie] +
			(_dataGhostState)[firstVelocity + 2*ie + 1]*(_dataGhostState)[firstVelocity + 2*ie + 1];
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

