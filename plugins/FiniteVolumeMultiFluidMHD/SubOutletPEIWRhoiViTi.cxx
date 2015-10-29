#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/SubOutletPEIWRhoiViTi.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
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

MethodCommandProvider<SubOutletPEIWRhoiViTi, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> subOutletPEIWRhoiViTiFVMCCProvider("SubOutletPEIWRhoiViTiFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void SubOutletPEIWRhoiViTi::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<CFreal> >("Pi","static pressure of the species");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletPEIWRhoiViTi::SubOutletPEIWRhoiViTi
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

SubOutletPEIWRhoiViTi::~SubOutletPEIWRhoiViTi()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletPEIWRhoiViTi::setup()
{
  FVMCC_BC::setup();
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
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

void SubOutletPEIWRhoiViTi::setGhostState(GeometricEntity *const face)
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
  CFreal nz = 0;
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;

  cf_assert(_updateVarSet.isNotNull());
  
  const CFreal bn = (_dataInnerState)[0]*nx + (_dataInnerState)[1]*ny;
  const CFreal en = (_dataInnerState)[3]*nx + (_dataInnerState)[4]*ny;  
//  const CFreal chi = _varSet->getModel()->getDivECleaningConst();

  (_dataGhostState)[0] = (_dataInnerState)[0] /*+ 2*bn*nx*/;	//Bx
  (_dataGhostState)[1] = (_dataInnerState)[1] /*+ 2*bn*ny*/;	//By
  (_dataGhostState)[2] = (_dataInnerState)[2] /*+ 2*bn*nz*/;	//Bz
  (_dataGhostState)[3] = (_dataInnerState)[3] /*- 2*en*nx*/;	//Ex
  (_dataGhostState)[4] = (_dataInnerState)[4] /*- 2*en*ny*/;	//Ey
  (_dataGhostState)[5] = (_dataInnerState)[5] /*- 2*en*nz*/;	//Ez
  (_dataGhostState)[6] = (_dataInnerState)[6];			//Psi
  (_dataGhostState)[7] = (_dataInnerState)[7];			//Phi
  
///MultiFluidMHD Subsonic Outlet Condition imposing pressure in 2D
 // here a fix is needed in order to have always ghostT > 0
 // if ghostT < 0  then the inner value is set
  const CFuint endEM = 8;
  const CFuint firstSpecies = _updateVarSet->getModel()->getFirstScalarVar(0);  
  const CFuint firstVelocity = _updateVarSet->getModel()->getFirstScalarVar(1);   
  const CFuint firstTemperature = _updateVarSet->getModel()->getFirstScalarVar(2); 
  const CFreal gamma = _updateVarSet->getModel()->getGamma();
  const CFreal Boltz = _updateVarSet->getModel()->getK();
  const CFreal molecularMass1 = _updateVarSet->getModel()->getMolecularMass1();
  const CFreal molecularMass2 = _updateVarSet->getModel()->getMolecularMass2();
  const CFreal molecularMass3 = _updateVarSet->getModel()->getMolecularMass3();
  const CFreal R1 = Boltz/molecularMass1;
  const CFreal R2 = Boltz/molecularMass2;
  const CFreal R3 = Boltz/molecularMass3;
  const CFreal gammaDivGammaMinus1 = gamma/(gamma - 1.0);
 
  (_dataGhostState)[endEM] = (_dataInnerState)[endEM]; 							//RHO
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    const CFreal Pi = _Pi[ie];
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
  
  // Temporary fix to one fluid 
  (_dataGhostState)[endEM] = (_dataGhostState)[firstTemperature + 1]/(R1*(_dataGhostState)[firstTemperature]);                //RHO is computed with the pressure
  (_dataGhostState)[firstSpecies] = 1;											      //y = 1
  const CFreal V2 = (_dataGhostState)[firstVelocity]*(_dataGhostState)[firstVelocity] +
                      (_dataGhostState)[firstVelocity + 1]*(_dataGhostState)[firstVelocity + 1];			      //Vx and Vy are correct;
  const CFreal rho =(_dataGhostState)[endEM];										      // we use this rho
  (_dataGhostState)[firstTemperature + 2] = sqrt(gamma*(_dataGhostState)[firstTemperature + 1]/rho);                         //ai
  (_dataGhostState)[firstTemperature + 3] = (0.5*rho*V2 +
                                            gammaDivGammaMinus1*(_dataGhostState)[firstTemperature + 1])/rho;                //Hi

  _updateVarSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
