#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/BCSubInletUVTPCWRhoiViTi.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Framework/PhysicalConsts.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/CFLog.hh"
#include "Common/NoSuchValueException.hh"
#include "FluxReconstructionMethod/StdSourceTerm.hh"

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

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<BCSubInletUVTPCWRhoiViTi, FluxReconstructionSolverData, BCStateComputer, FluxReconstructionMultiFluidMHDModule> BCsubInletUVTPCWRhoiViTiProvider("BCSubInletUVTPCWRhoiViTi");

//////////////////////////////////////////////////////////////////////////////

void BCSubInletUVTPCWRhoiViTi::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");      
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletUVTPCWRhoiViTi::BCSubInletUVTPCWRhoiViTi(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_uvT(),
  m_useFunction(false),
  m_bCoord(),
  m_vars()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

   m_functions = std::vector<std::string>();
   setParameter("Def",&m_functions);

   m_vars = std::vector<std::string>();
   setParameter("Vars",&m_vars);   
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletUVTPCWRhoiViTi::~BCSubInletUVTPCWRhoiViTi()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletUVTPCWRhoiViTi::setup()
{
  BCStateComputer::setup();
  
//   cout << "SubInletUVTPCWRhoiViTi::setup" << endl; 
  
  m_needsSpatCoord = true;
  
  
  m_bCoord.resize(PhysicalModelStack::getActive()->getDim());
  
//   _vars.resize(PhysicalModelStack::getActive()->getDim());
  
  m_varSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  
  // resize the physical data for internal and ghost solution points
  m_varSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_varSet->getModel()->resizePhysicalData(m_intSolPhysData);

  const CFuint nbSpecies = m_varSet->getModel()->getNbScalarVars(0); 
  
  m_uvT.resize(3*nbSpecies);  
 
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletUVTPCWRhoiViTi::configure ( Config::ConfigArgs& args )
{
  using namespace COOLFluiD::Framework;

  BCStateComputer::configure(args);

    m_vFunction.setFunctions(m_functions);
    m_vFunction.setVariables(m_vars);
    try {
      m_vFunction.parse();
      m_useFunction = true;
    }
    catch (Common::ParserException& e) {
      CFout << e.what() << "\n";
      throw; // retrow the exception to signal the error to the user
    }

}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletUVTPCWRhoiViTi::computeGhostStates(const vector< State* >& intStates,
                                                    vector< State* >& ghostStates,
                                                    const std::vector< RealVector >& normals,
                                                    const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbSpecies = m_varSet->getModel()->getNbScalarVars(0);
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

//  CFLog(VERBOSE, "\n\n\n\n\n  nbSpecies Size = " << nbSpecies << "\n\n\n\n\n");
//  CFLog(VERBOSE, "\n\n\n\n\n  nbrStates Size = " << nbrStates << "\n\n\n\n\n");

//  CFLog(VERBOSE, "\n\n\n Error 1 \n");

  // get some physical data from the model
  const CFreal gamma = m_varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);


  ///MultiFluidMHD Subsonic Outlet Condition imposing pressure in 2D
  // here a fix is needed in order to have always ghostT > 0
  // if ghostT < 0  then the inner value is set
  const CFuint endEM = 8;
  
  //CFLog(VERBOSE, "\n\n\n\n\n  m_ghostSolPhysData Size = " << m_ghostSolPhysData.size() << "\n\n\n\n\n");
  //CFLog(VERBOSE, "\n\n\n\n\n  m_intSolPhysData Size = " << m_intSolPhysData.size() << "\n\n\n\n\n");
  
  //const CFuint tgAlpha = tan(m_alpha);
  

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
      // dereference states
      State& intState   = (*intStates[iState]);
      State& ghostState = (*ghostStates[iState]);

//  CFLog(VERBOSE, "\n\n\n Error 2 \n");

      // set the physical data starting from the inner state
//      m_varSet->computePhysicalData(intState, m_intSolPhysData);      

      //m_varSet->computePhysicalData(*((*m_cellStates)[iState]), m_intSolPhysData);

    //const CFuint stateID = (*intState)[iState]->getLocalID();

      // normal
      const RealVector& normal = normals[iState];

//  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  CFreal nx = normal[XX];
  CFreal ny = normal[YY];
  CFreal nz = 0;
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;

  cf_assert(m_varSet.isNotNull());
  
  const CFreal bn = (intState)[0]*nx + (intState)[1]*ny;
  const CFreal en = (intState)[3]*nx + (intState)[4]*ny;  
//  const CFreal chi = _varSet->getModel()->getDivECleaningConst();

  (ghostState)[0] = (intState)[0] - 2*bn*nx;	//Bx
  (ghostState)[1] = (intState)[1] - 2*bn*ny;	//By
  (ghostState)[2] = (intState)[2] - 2*bn*nz;	//Bz
  (ghostState)[3] = -(intState)[3] + 2*en*nx;	//Ex
  (ghostState)[4] = -(intState)[4] + 2*en*ny;	//Ey
  (ghostState)[5] = -(intState)[5] + 2*en*nz;	//Ez
  (ghostState)[6] = (intState)[6];			//Psi
  (ghostState)[7] = -(intState)[7];			//Phi

//   std::cout << "SubInletUVTPCWRhoiViTi::setGhostState after Maxwell" <<"\n"; 

///MultiFluidMHD No Slip Isothermal Condition in 2D
 // here a fix is needed in order to have always ghostT > 0
 // if ghostT < 0  then the inner value is set
 const CFuint endEM = 8;


 //set the densities
 for (CFuint i = 0 ; i < nbSpecies; i++){
   (ghostState)[endEM + i] = (intState)[endEM + i];
 }
 
  //CFLog(VERBOSE, "\n\n\n Error 3 \n");

 if(m_useFunction){
    // coordinate of the boundary point
   //m_bCoord = intState.getCoordinates();
               //ghostState.getCoordinates());

//  CFLog(VERBOSE, "\n\n\n\n\n intState.getCoordinates() Value = " << intState.getCoordinates() << "\n\n\n\n\n");
//  CFLog(VERBOSE, "\n\n\n\n\n ghostState.getCoordinates() Value = " << ghostState.getCoordinates() << "\n\n\n\n\n");

//   m_bCoord = ((*m_cellStates)[nbSpecies]->getCoordinates()); 
//               // + ghostState.getCoordinates());
//   m_bCoord *= 0.5;

   // (*ghostState) = 2*bcState - (*innerState)
   m_vFunction.evaluate(coords[iState], m_uvT);
   //CFLog(INFO, "coords: " << coords[iState] << ", m_uvT: " << m_uvT << "\n");

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
   const CFreal Ui = m_uvT[i];
   const CFreal Vi = m_uvT[nbSpecies + i];
//    std::cout << "SubInletUVTPCWRhoiViTi::setGhostState Ui = "<< Ui <<"\n"; 
//    std::cout << "SubInletUVTPCWRhoiViTi::setGhostState Vi = "<< Vi <<"\n"; 

   (ghostState)[endEM + nbSpecies + 2*i] = 2*Ui -(intState)[endEM + nbSpecies + 2*i];
   //CFLog(INFO, "Ui: " << Ui << ", ghostU: " << (ghostState)[endEM + nbSpecies + 2*i] << "\n");
   (ghostState)[endEM + nbSpecies + 2*i + 1] = 2*Vi  -(intState)[endEM + nbSpecies + 2*i + 1];   
 } 
//   std::cout << "SubInletUVTPCWRhoiViTi::setGhostState after Velocities" <<"\n"; 

 //set the Temperatures
 for (CFuint i = 0 ; i < nbSpecies; i++){
   const CFreal Ti = m_uvT[2*nbSpecies + i];   
   (ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = 2.*Ti - (intState)[endEM + nbSpecies + 2*nbSpecies + i];
   cf_assert(2.*Ti - (intState)[endEM + nbSpecies + 2*nbSpecies + i] > 0.); 
 }  

    //m_varSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);

   //CFLog(INFO, "ghost inside inlet: " << ghostState << "\n");
    //CFLog(VERBOSE, "\n\n\n\n\n Physical Data Value = " << m_ghostSolPhysData << "\n\n\n\n\n");
 }
  
}
//////////////////////////////////////////////////////////////////////////////

void BCSubInletUVTPCWRhoiViTi::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                            std::vector< std::vector< RealVector* > >& ghostGrads,
                                            const std::vector< RealVector >& normals,
                                            const std::vector< RealVector >& coords)
{

  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // number of gradient variables
  cf_assert(nbrStateGrads > 0);
  const CFuint nbrGradVars = intGrads[0].size();

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
    {
      *ghostGrads[iState][iGradVar] = *intGrads[iState][iGradVar];
    }
  }

/*
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // number of gradient variables
  cf_assert(nbrStateGrads > 0);

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    cf_assert(intGrads[iState].size() == 4);

    // normal
    const RealVector& normal = normals[iState];

    // pressure
    RealVector& presGradI = *intGrads  [iState][0];
    RealVector& presGradG = *ghostGrads[iState][0];
    const CFreal nPresGrad = (presGradI[XX]*normal[XX] + presGradI[YY]*normal[YY]);
    presGradG = presGradI - 2.0*nPresGrad*normal;

    // temperature
    RealVector& tempGradI = *intGrads  [iState][4];
    RealVector& tempGradG = *ghostGrads[iState][4];
    const CFreal nTempGrad = (tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY]);
    tempGradG = tempGradI - 2.0*nTempGrad*normal;

    // velocity
    RealVector& uGradI = *intGrads  [iState][1];
    RealVector& uGradG = *ghostGrads[iState][1];
    RealVector& vGradI = *intGrads  [iState][2];
    RealVector& vGradG = *ghostGrads[iState][2];
    //RealVector& wGradI = *intGrads  [iState][3];
    //RealVector& wGradG = *ghostGrads[iState][3];

    // internal normal and tangential component
    const RealVector velocNGradI = uGradI*normal[XX] + vGradI*normal[YY];
    const RealVector uTGradI = uGradI - velocNGradI*normal[XX];
    const RealVector vTGradI = vGradI - velocNGradI*normal[YY];
    //const RealVector wTGradI = wGradI - velocNGradI*normal[ZZ];

    // ghost normal and tangential component
    const RealVector velocNGradG = velocNGradI;
    RealVector velocTGradNI(3);
    velocTGradNI[XX] = uTGradI[XX]*normal[XX] + uTGradI[YY]*normal[YY]; // + uTGradI[ZZ]*normal[ZZ];
    velocTGradNI[YY] = vTGradI[XX]*normal[XX] + vTGradI[YY]*normal[YY]; // + vTGradI[ZZ]*normal[ZZ];
    //velocTGradNI[ZZ] = wTGradI[XX]*normal[XX] + wTGradI[YY]*normal[YY] + wTGradI[ZZ]*normal[ZZ];
    const RealVector uTGradG = uTGradI - 2.0*velocTGradNI[XX]*normal;
    const RealVector vTGradG = vTGradI - 2.0*velocTGradNI[YY]*normal;
    //const RealVector wTGradG = wTGradI - 2.0*velocTGradNI[ZZ]*normal;

    // compute ghost velocity gradients
    uGradG = uTGradG + velocNGradG*normal[XX];
    vGradG = vTGradG + velocNGradG*normal[YY];
    //wGradG = wTGradG + velocNGradG*normal[ZZ];
  }
*/
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

