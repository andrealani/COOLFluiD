#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/SubInletUVTEIWRhoiViTi.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
// #include "Framework/MeshData.hh"
  
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace COOLFluiD::Physics::Maxwell;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubInletUVTEIWRhoiViTi, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> subInletUVTEIWRhoiViTiFVMCCProvider("SubInletUVTEIWRhoiViTiFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void SubInletUVTEIWRhoiViTi::defineConfigOptions(Config::OptionList& options)
{
//    options.addConfigOption< std::vector<CFreal> >("Vx","x velocity of the species");
//    options.addConfigOption< std::vector<CFreal> >("Vy","y velocity of the species");
//    options.addConfigOption< std::vector<CFreal> >("T","static temperature of the species");
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");      
}

//////////////////////////////////////////////////////////////////////////////

SubInletUVTEIWRhoiViTi::SubInletUVTEIWRhoiViTi
(const std::string& name) :
  FVMCC_BC(name),
  _uvT(),
  _useFunction(false),
  _updateVarSet(CFNULL),
  _bCoord(),
  _vars()
{
   addConfigOptionsTo(this);
   
//    _ui = std::vector<CFreal>();
//    setParameter("Vx",&_ui);
//    
//    _vi = std::vector<CFreal>();
//    setParameter("Vy",&_vi);
//    
//    _Ti = std::vector<CFreal>();
//    setParameter("Ti",&_Ti);
   
   _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

   _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);   
}
      
//////////////////////////////////////////////////////////////////////////////

SubInletUVTEIWRhoiViTi::~SubInletUVTEIWRhoiViTi()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletUVTEIWRhoiViTi::setup()
{
  FVMCC_BC::setup();
  
//   cout << "SubInletUVTEIWRhoiViTi::setup" << endl; 
  
  
  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
  
//   _vars.resize(PhysicalModelStack::getActive()->getDim());
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0); 
  
  _uvT.resize(3*nbSpecies);  
 
}

//////////////////////////////////////////////////////////////////////////////

void SubInletUVTEIWRhoiViTi::configure ( Config::ConfigArgs& args )
{
  using namespace COOLFluiD::Framework;

  FVMCC_BC::configure(args);

//   if(!_functions.empty())
//     {
    _vFunction.setFunctions(_functions);
    _vFunction.setVariables(_vars);
    try {
      _vFunction.parse();
      _useFunction = true;
    }
    catch (Common::ParserException& e) {
      CFout << e.what() << "\n";
      throw; // retrow the exception to signal the error to the user
    }
//   }
}

//////////////////////////////////////////////////////////////////////////////

void SubInletUVTEIWRhoiViTi::setGhostState(GeometricEntity *const face)
{
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);
  
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  
  ///Maxwell Equations Perfectly Conducting Wall Condition
//   std::cout << "SubInletUVTEIWRhoiViTi::setGhostState before assignment" <<"\n";
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
  
  const CFreal bn = (*innerState)[0]*nx + (*innerState)[1]*ny;
  const CFreal en = (*innerState)[3]*nx + (*innerState)[4]*ny;  
//  const CFreal chi = _varSet->getModel()->getDivECleaningConst();

  (*ghostState)[0] = (*innerState)[0] /*+ 2*bn*nx*/;	//Bx
  (*ghostState)[1] = (*innerState)[1] /*+ 2*bn*ny*/;	//By
  (*ghostState)[2] = (*innerState)[2] /*+ 2*bn*nz*/;	//Bz
  (*ghostState)[3] = (*innerState)[3] /*- 2*en*nx*/;	//Ex
  (*ghostState)[4] = (*innerState)[4] /*- 2*en*ny*/;	//Ey
  (*ghostState)[5] = (*innerState)[5] /*- 2*en*nz*/;	//Ez
  (*ghostState)[6] = (*innerState)[6];			//Psi
  (*ghostState)[7] = (*innerState)[7];			//Phi

//   std::cout << "SubInletUVTEIWRhoiViTi::setGhostState after Maxwell" <<"\n"; 

///MultiFluidMHD No Slip Isothermal Condition in 2D
 // here a fix is needed in order to have always ghostT > 0
 // if ghostT < 0  then the inner value is set
 const CFuint endEM = 8;
 
 //set the densities
 for (CFuint i = 0 ; i < nbSpecies; i++){
   (*ghostState)[endEM + i] = (*innerState)[endEM + i];
 }
 
 if(_useFunction){
    // coordinate of the boundary point
   _bCoord = (innerState->getCoordinates() +
               ghostState->getCoordinates());
   _bCoord *= 0.5;

   // (*ghostState) = 2*bcState - (*innerState)
   _vFunction.evaluate(_bCoord, _uvT);

//    _uvT[0] /= _uvTRef[0];
//    _uvT[1] /= _uvTRef[1];
//    _uvT[2] /= _uvTRef[2];
 }
//  else {
//    for (CFuint i = 0 ; i < nbSpecies; i++){
// 	 const CFreal Ui = _ui[i];
// 	 const CFreal Vi = _vi[i];
// 	 const CFreal Ti = _Ti[i];
//      _uvT[i] = Ui;
//      _uvT[nbSpecies + i] = Vi;
//      _uvT[2*nbSpecies + i] = Ti;
//    }
//  }
 
 //set the Velocities
 for (CFuint i = 0 ; i < nbSpecies; i++){
   const CFreal Ui = _uvT[i];
   const CFreal Vi = _uvT[nbSpecies + i];
//    std::cout << "SubInletUVTEIWRhoiViTi::setGhostState Ui = "<< Ui <<"\n"; 
//    std::cout << "SubInletUVTEIWRhoiViTi::setGhostState =>     Vi        = "<< Vi <<"\n"; 

   (*ghostState)[endEM + nbSpecies + 2*i] = 2.*Ui -(*innerState)[endEM + nbSpecies + 2*i];
   (*ghostState)[endEM + nbSpecies + 2*i + 1] = 2.*Vi  -(*innerState)[endEM + nbSpecies + 2*i + 1];
//    cf_assert(Vi == 0.);
   
//    cf_assert((*ghostState)[endEM + nbSpecies + 2*i + 1] + (*innerState)[endEM + nbSpecies + 2*i + 1] != 0.);
   
//    std::cout << "SubInletUVTEIWRhoiViTi::setGhostState => (*ghostState) = "<<    (*ghostState)[endEM + nbSpecies + 2*i + 1] <<"\n"; 
//    std::cout << "SubInletUVTEIWRhoiViTi::setGhostState => (*innerState) = "<<    (*innerState)[endEM + nbSpecies + 2*i + 1] <<"\n"; 
   
 } 
//   std::cout << "SubInletUVTEIWRhoiViTi::setGhostState after Velocities" <<"\n"; 
 //set the Temperatures
 for (CFuint i = 0 ; i < nbSpecies; i++){
   const CFreal Ti = _uvT[2*nbSpecies + i];   
   (*ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = 2.*Ti - (*innerState)[endEM + nbSpecies + 2*nbSpecies + i];
   cf_assert(2.*Ti - (*innerState)[endEM + nbSpecies + 2*nbSpecies + i] > 0.); 
 }  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
