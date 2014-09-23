#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/NoSlipWallIsothermalEIWRhoiViTi.hh"
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

MethodCommandProvider<NoSlipWallIsothermalEIWRhoiViTi, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> noSlipWallIsothermalEIWRhoiViTiFVMCCProvider("NoSlipWallIsothermalEIWRhoiViTiFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalEIWRhoiViTi::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<CFreal> >("TWall","Wall temperatures of the species");
  options.addConfigOption< CFreal > ("Ez0","Imposed Ez0");
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermalEIWRhoiViTi::NoSlipWallIsothermalEIWRhoiViTi
(const std::string& name) :
  FVMCC_BC(name),
  _updateVarSet(CFNULL)
{
  addConfigOptionsTo(this);
  _wallTemp = std::vector<CFreal>();
  setParameter("TWall",&_wallTemp);
   
  _Ez0 = -0.1;
  setParameter("Ez0",&_Ez0);
}
      
//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermalEIWRhoiViTi::~NoSlipWallIsothermalEIWRhoiViTi()
{
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalEIWRhoiViTi::setup()
{
  FVMCC_BC::setup();
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  
  
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);  
  for (CFuint i = 0; i < nbSpecies; ++i) {
    cf_assert(_wallTemp[i] >= 0.0);
  }

  // adimensionalize the temperature
//   _wallTemp /= _updateVarSet->getModel()->getTempRef();
}
      
//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalEIWRhoiViTi::setGhostState(GeometricEntity *const face)
{
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);
  
  State *const innerState = face->getState(0); 
  State *const ghostState = face->getState(1);
  
  
  ///Maxwell Equations Perfectly Conducting Wall Condition
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
  const CFreal muZero = 1.2566370614359e-6;
  const CFreal omega = 1.82e5; //14857.14286;
  const CFreal Ez0 = _Ez0;  
   
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = 0;
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;
  
//   const CFreal dy = ghostState->getCoordinates()[YY] - innerState->getCoordinates()[YY];
  const CFreal dr = MathFunctions::getDistance(ghostState->getCoordinates(),
						  innerState->getCoordinates());
  //cout<<"dr = "<< dr <<"\n";

  cf_assert(_updateVarSet.isNotNull());
  
  const CFreal bn = (*innerState)[0]*nx + (*innerState)[1]*ny;
  const CFreal en = (*innerState)[3]*nx + (*innerState)[4]*ny;  
//  const CFreal chi = _varSet->getModel()->getDivECleaningConst();

  (*ghostState)[0] = (*innerState)[0] - dr*muZero*omega*Ez0;	//Bx
  (*ghostState)[1] = -(*innerState)[1]; //By
  (*ghostState)[2] = (*innerState)[2]; //Bz
  (*ghostState)[3] = -(*innerState)[3]; //Ex
  (*ghostState)[4] = -(*innerState)[4]; //Ey
  (*ghostState)[5] = -(*innerState)[5]; //Ez
  (*ghostState)[6] = -(*innerState)[6];	//Psi
  (*ghostState)[7] = -(*innerState)[7];	//Phi
  
///MultiFluidMHD No Slip Isothermal Condition in 2D
 // here a fix is needed in order to have always ghostT > 0
 // if ghostT < 0  then the inner value is set
 const CFuint endEM = 8;
 
 //set the densities
 for (CFuint i = 0 ; i < nbSpecies; i++){
   (*ghostState)[endEM + i] = (*innerState)[endEM + i];
 }
 
 //set the Velocities
 for (CFuint i = 0 ; i < nbSpecies; i++){
   (*ghostState)[endEM + nbSpecies + 2*i] = -(*innerState)[endEM + nbSpecies + 2*i];
   (*ghostState)[endEM + nbSpecies + 2*i + 1] = -(*innerState)[endEM + nbSpecies + 2*i + 1];   
 } 
 
 //set the Temperatures
 for (CFuint i = 0 ; i < nbSpecies; i++){
   const CFreal Twall = _wallTemp[i];
   (*ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = 2.*Twall - (*innerState)[endEM + nbSpecies + 2*nbSpecies + i];
   cf_assert(2.*Twall - (*innerState)[endEM + nbSpecies + 2*nbSpecies + i] > 0.); 
 } 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
