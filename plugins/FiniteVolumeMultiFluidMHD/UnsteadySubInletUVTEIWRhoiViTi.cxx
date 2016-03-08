#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/UnsteadySubInletUVTEIWRhoiViTi.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalConsts.hh"
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

MethodCommandProvider<UnsteadySubInletUVTEIWRhoiViTi, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> UnsteadySubInletUVTEIWRhoiViTiFVMCCProvider("UnsteadySubInletUVTEIWRhoiViTiFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletUVTEIWRhoiViTi::defineConfigOptions(Config::OptionList& options)
{
//    options.addConfigOption< std::vector<CFreal> >("Vx","x velocity of the species");
//    options.addConfigOption< std::vector<CFreal> >("Vy","y velocity of the species");
//    options.addConfigOption< std::vector<CFreal> >("T","static temperature of the species");
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");      
}

//////////////////////////////////////////////////////////////////////////////

UnsteadySubInletUVTEIWRhoiViTi::UnsteadySubInletUVTEIWRhoiViTi
(const std::string& name) :
  FVMCC_BC(name),
  _uvT(),
  _useFunction(false),
  _updateVarSet(CFNULL),
  _bCoord()
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

UnsteadySubInletUVTEIWRhoiViTi::~UnsteadySubInletUVTEIWRhoiViTi()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletUVTEIWRhoiViTi::setup()
{
  FVMCC_BC::setup();
  
//   cout << "UnsteadySubInletUVTEIWRhoiViTi::setup" << endl;
  
  
  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
  
  _variables.resize(PhysicalModelStack::getActive()->getDim() + 1);
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();

  _updateVarSet->getModel()->resizePhysicalData(_physicalData);  

  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0); 
  
  _uvT.resize(3*nbSpecies + 3);  
 
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletUVTEIWRhoiViTi::configure ( Config::ConfigArgs& args )
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

void UnsteadySubInletUVTEIWRhoiViTi::setGhostState(GeometricEntity *const face)
{
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);
  
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  // compute the physical data
  _updateVarSet->computePhysicalData(*innerState, _physicalData); 
 
  ///Maxwell Equations Perfectly Conducting Wall Condition
//   std::cout << "UnsteadySubInletUVTEIWRhoiViTi::setGhostState before assignment" <<"\n";
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

 if(_useFunction){
    // coordinate of the boundary point
    _bCoord = (innerState->getCoordinates() +
                 ghostState->getCoordinates());
    _bCoord *= 0.5;
    
    const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
    for(CFuint iDim = 0; iDim < nbDim; iDim++)
    {
      _variables[iDim] = _bCoord[iDim];
    }
    _variables[PhysicalModelStack::getActive()->getDim()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();
    _vFunction.evaluate(_variables, _uvT);
    
  }
  
  const CFreal BxBound = _uvT[3*nbSpecies];
  const CFreal ByBound = _uvT[3*nbSpecies + 1];
  const CFreal EzBound = _uvT[3*nbSpecies + 2];
  (*ghostState)[0] = (*innerState)[0] /*+ 2*bn*nx*/;	//Bx is harcoded to be extrapolated from inside
  (*ghostState)[1] = 2*ByBound - (*innerState)[1] /*+ 2*bn*ny*/;	//By
  (*ghostState)[2] = (*innerState)[2] /*+ 2*bn*nz*/;	//Bz
  (*ghostState)[3] = (*innerState)[3] /*- 2*en*nx*/;	//Ex
  (*ghostState)[4] = (*innerState)[4] /*- 2*en*ny*/;	//Ey
  (*ghostState)[5] = (*innerState)[5] /*- 2*en*nz*/;	//Ez is harcoded to be extrapolated from inside
  (*ghostState)[6] = -(*innerState)[6];			//Psi
  (*ghostState)[7] = (*innerState)[7];			//Phi

//   std::cout << "UnsteadySubInletUVTEIWRhoiViTi::setGhostState after Maxwell" <<"\n";

///MultiFluidMHD No Slip Isothermal Condition in 2D
 // here a fix is needed in order to have always ghostT > 0
 // if ghostT < 0  then the inner value is set
 const CFuint endEM = 8;
 
 //previous set for the densities - gives wrong pressure at the boundary
// for (CFuint i = 0 ; i < nbSpecies; i++){
//   (*ghostState)[endEM + i] = (*innerState)[endEM + i];
// }

 const CFreal m_e = _updateVarSet->getModel()->getMolecularMass1();
 const CFreal m_n = _updateVarSet->getModel()->getMolecularMass2();
 const CFreal m_p = _updateVarSet->getModel()->getMolecularMass3();
 
 // AL: gory fix here: assume 2-fluid with protons and neutrals
 // This needs to be generalized by AAL  
 RealVector mi(nbSpecies); 
 if (nbSpecies == 2 && _updateVarSet->getModel()->isLeake()) {
   mi[0] = m_p; mi[1] = m_n;
 }
 else {
   CFLog(ERROR, "UnsteadySubInletUVTEIWRhoiViTi::setGhostState() => only implemented for Leake model!");
   abort();
 }
 
 // new: set the temperatures and densities
// const CFreal rhoInner = _physicalData[EulerMFMHDTerm::RHO];
 const CFuint firstTemperature = _updateVarSet->getModel()->getFirstScalarVar(2);
 for (CFuint i = 0 ; i < nbSpecies; i++){
   const CFreal Ti = _uvT[2*nbSpecies + i];
//   const CFreal rhoiGhost = _uvT[i];
   const CFreal Tighost = 2.*Ti - (*innerState)[endEM + 3*nbSpecies + i];
   (*ghostState)[endEM + 3*nbSpecies + i] = Tighost;
   cf_assert(2.*Ti - (*innerState)[endEM + 3*nbSpecies + i] > 0.); 
   
   //const CFreal yiInner  = (*innerState)[endEM + i]/rhoInner;
   //const CFreal yiGhost  = 2.*_uvT[i] - yiInner;
   //we want p_i ghost = p_i inner, so we set the density
   const CFreal piGhost  = _physicalData[firstTemperature + i*4 + 1];  // imposing Pghost = Pinner
   const CFreal R_gas = PhysicalConsts::Boltzmann()/mi[i];
   //const CFreal rhoGhost = piGhost/(yiGhost*R_gas*Tighost);
   if(i == 0){(*ghostState)[endEM + i] = piGhost/(2*R_gas*Tighost);}
   else{   (*ghostState)[endEM + i] = piGhost/(R_gas*Tighost);} // imposing rho_g = p_g/RT_g
//   std::cout << "Ti = " << Ti <<"\n";
//   std::cout << "Tighost = " << Tighost <<"\n";
//   std::cout << "rhoiGhost = " << (*ghostState)[endEM + i] <<"\n";
//   std::cout << "rhoiInner = " << (*innerState)[endEM + i] <<"\n";
//   std::cout << "piGhost = " << piGhost << "\n";
}
//   abort(); 


 //set the Velocities
 for (CFuint i = 0 ; i < nbSpecies; i++){
   const CFreal Ui = _uvT[i];
   const CFreal Vi = _uvT[nbSpecies + i];
//    std::cout << "UnsteadySubInletUVTEIWRhoiViTi::setGhostState Ui = "<< Ui <<"\n";
//    std::cout << "UnsteadySubInletUVTEIWRhoiViTi::setGhostState =>     Vi        = "<< Vi <<"\n";

   (*ghostState)[endEM + nbSpecies + 2*i] = 2.*Ui -(*innerState)[endEM + nbSpecies + 2*i];
   (*ghostState)[endEM + nbSpecies + 2*i + 1] = 2.*Vi  -(*innerState)[endEM + nbSpecies + 2*i + 1];

//    cf_assert(Vi == 0.);
   
//    cf_assert((*ghostState)[endEM + nbSpecies + 2*i + 1] + (*innerState)[endEM + nbSpecies + 2*i + 1] != 0.);
   
    //std::cout << "UnsteadySubInletUVTEIWRhoiViTi::setGhostState => (*ghostState) = "<<    (*ghostState)[endEM + nbSpecies + 2*i + 1] <<"\n";
    //std::cout << "UnsteadySubInletUVTEIWRhoiViTi::setGhostState => (*innerState) = "<<    (*innerState)[endEM + nbSpecies + 2*i + 1] <<"\n";
   
 } 
//   std::cout << "UnsteadySubInletUVTEIWRhoiViTi::setGhostState after Velocities" <<"\n";
 //set the Temperatures
 //for (CFuint i = 0 ; i < nbSpecies; i++){
 //  const CFreal Ti = _uvT[2*nbSpecies + i];
 //  (*ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = 2.*Ti - (*innerState)[endEM + nbSpecies + 2*nbSpecies + i];
 //  cf_assert(2.*Ti - (*innerState)[endEM + nbSpecies + 2*nbSpecies + i] > 0.); 
 //}  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
