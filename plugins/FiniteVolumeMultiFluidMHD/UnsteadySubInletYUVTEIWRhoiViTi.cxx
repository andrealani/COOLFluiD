#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/UnsteadySubInletYUVTEIWRhoiViTi.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalConsts.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
  
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

MethodCommandProvider<UnsteadySubInletYUVTEIWRhoiViTi, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> 
UnsteadySubInletYUVTEIWRhoiViTiFVMCCProvider("UnsteadySubInletYUVTEIWRhoiViTiFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletYUVTEIWRhoiViTi::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");      
}

//////////////////////////////////////////////////////////////////////////////

UnsteadySubInletYUVTEIWRhoiViTi::UnsteadySubInletYUVTEIWRhoiViTi
(const std::string& name) :
  FVMCC_BC(name),
  _yuvT(),
  _useFunction(false),
  _updateVarSet(CFNULL),
  _bCoord()
{
  addConfigOptionsTo(this);
  
  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);
  
  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);   
}
      
//////////////////////////////////////////////////////////////////////////////

UnsteadySubInletYUVTEIWRhoiViTi::~UnsteadySubInletYUVTEIWRhoiViTi()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletYUVTEIWRhoiViTi::setup()
{
  FVMCC_BC::setup();
  
  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
  
  _variables.resize(PhysicalModelStack::getActive()->getDim() + 1);

  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  
  _updateVarSet->getModel()->resizePhysicalData(_physicalData);
  
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0); 
  // for each species we will read mass fraction, velocity components and temperature
  _yuvT.resize(4*nbSpecies + 3);  
}
      
//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletYUVTEIWRhoiViTi::configure ( Config::ConfigArgs& args )
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

void UnsteadySubInletYUVTEIWRhoiViTi::setGhostState(GeometricEntity *const face)
{
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);
  
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  // compute the physical data
  _updateVarSet->computePhysicalData(*innerState, _physicalData);
  
  ///Maxwell Equations Perfectly Conducting Wall Condition
  //   std::cout << "UnsteadySubInletYUVTEIWRhoiViTi::setGhostState before assignment" <<"\n";
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
    _bCoord = (innerState->getCoordinates() + ghostState->getCoordinates());
    _bCoord *= 0.5;
    
    const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
    for(CFuint iDim = 0; iDim < nbDim; iDim++) {
      _variables[iDim] = _bCoord[iDim];
    }
    _variables[PhysicalModelStack::getActive()->getDim()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();
    _vFunction.evaluate(_variables, _yuvT);
  }
  
  const CFreal BxBound = _yuvT[4*nbSpecies];
  const CFreal ByBound = _yuvT[4*nbSpecies + 1];
  const CFreal EzBound = _yuvT[4*nbSpecies + 2];
  (*ghostState)[0] = (*innerState)[0] /*+ 2*bn*nx*/;	//Bx is harcoded to be extrapolated from inside
  (*ghostState)[1] = 2*ByBound - (*innerState)[1] /*+ 2*bn*ny*/;	//By
  (*ghostState)[2] = (*innerState)[2] /*+ 2*bn*nz*/;	//Bz
  (*ghostState)[3] = (*innerState)[3] /*- 2*en*nx*/;	//Ex
  (*ghostState)[4] = (*innerState)[4] /*- 2*en*ny*/;	//Ey
  (*ghostState)[5] = (*innerState)[5] /*- 2*en*nz*/;	//Ez is harcoded to be extrapolated from inside
  (*ghostState)[6] = -(*innerState)[6];			//Psi
  (*ghostState)[7] = (*innerState)[7];			//Phi

//   std::cout << "UnsteadySubInletYUVTEIWRhoiViTi::setGhostState after Maxwell" <<"\n";

///MultiFluidMHD No Slip Isothermal Condition in 2D
 // here a fix is needed in order to have always ghostT > 0
 // if ghostT < 0  then the inner value is set
 const CFuint endEM = 8;
 
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
   CFLog(ERROR, "UnsteadySubInletYUVTEIWRhoiViTi::setGhostState() => only implemented for Leake model!");
   abort();
 }
 
 // set the temperatures and partial densities
 const CFreal rhoInner = _physicalData[EulerMFMHDTerm::RHO];
 const CFuint firstTemperature = _updateVarSet->getModel()->getFirstScalarVar(2);
 for (CFuint i = 0 ; i < nbSpecies; i++){
   const CFreal Ti = _yuvT[3*nbSpecies + i];
   const CFreal Tghost = 2.*Ti - (*innerState)[endEM + 3*nbSpecies + i];
   (*ghostState)[endEM + 3*nbSpecies + i] = Tghost;
   cf_assert(2.*Ti - (*innerState)[endEM + 3*nbSpecies + i] > 0.); 
   
   const CFreal yiInner  = (*innerState)[endEM + i]/rhoInner;
   const CFreal yiGhost  = 2.*_yuvT[i] - yiInner;
   // p_i ghost = p_i inner
   const CFreal piGhost  = _physicalData[firstTemperature + i*4 + 1];  // TO BE CHECKED!!!
   const CFreal R_gas = PhysicalConsts::Boltzmann()/mi[i];
   const CFreal rhoGhost = piGhost/(yiGhost*R_gas*Tghost);
   (*ghostState)[endEM + i] = rhoGhost*yiGhost;
 }
 
 //set the Velocities
 for (CFuint i = 0 ; i < nbSpecies; i++){
   const CFreal Ui = _yuvT[nbSpecies + i];
   const CFreal Vi = _yuvT[2*nbSpecies + i];
   //    std::cout << "UnsteadySubInletYUVTEIWRhoiViTi::setGhostState Ui = "<< Ui <<"\n";
   //    std::cout << "UnsteadySubInletYUVTEIWRhoiViTi::setGhostState =>     Vi        = "<< Vi <<"\n";
   (*ghostState)[endEM + nbSpecies + 2*i] = 2.*Ui -(*innerState)[endEM + nbSpecies + 2*i];
   (*ghostState)[endEM + nbSpecies + 2*i + 1] = 2.*Vi  -(*innerState)[endEM + nbSpecies + 2*i + 1];
   
   //    cf_assert(Vi == 0.);
   
   //    cf_assert((*ghostState)[endEM + nbSpecies + 2*i + 1] + (*innerState)[endEM + nbSpecies + 2*i + 1] != 0.);
   
   //std::cout << "UnsteadySubInletYUVTEIWRhoiViTi::setGhostState => (*ghostState) = "<<    (*ghostState)[endEM + nbSpecies + 2*i + 1] <<"\n";
   //std::cout << "UnsteadySubInletYUVTEIWRhoiViTi::setGhostState => (*innerState) = "<<    (*innerState)[endEM + nbSpecies + 2*i + 1] <<"\n";
 } 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
