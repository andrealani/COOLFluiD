#include "Framework/MethodStrategyProvider.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/BCNoSlipWallHeatFluxMF2D.hh"

#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace COOLFluiD::Physics::Maxwell;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCNoSlipWallHeatFluxMF2D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMultiFluidMHDModule >
  BCNoSlipWallHeatFluxNMF2DProvider("BCNoSlipWallHeatFluxMF2D");

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallHeatFluxMF2D::BCNoSlipWallHeatFluxMF2D(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;
  
  addConfigOptionsTo(this);

  _wallTemp = std::vector<CFreal>();
  setParameter("T",&_wallTemp);

//  m_wallT = 0.0; //std::vector<CFreal>(); 
//  setParameter("T",&m_wallT);

  m_wallQ = 0.0; //std::vector<CFreal>();
   setParameter("q",&m_wallQ);

  _Ez0 = -0.1;
  setParameter("Ez0",&_Ez0);

  m_heatFlux= true;
   setParameter("HeatFlux",&m_heatFlux);
   
  m_changeToIsoT = MathTools::MathConsts::CFuintMax();
   setParameter("ChangeToIsoT",&m_changeToIsoT);
}

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallHeatFluxMF2D::~BCNoSlipWallHeatFluxMF2D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxMF2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<CFreal> >("T","wall static temperature");
  options.addConfigOption< CFreal > ("Ez0","Imposed Ez0");
  options.addConfigOption< CFreal >("q","wall heat flux");
  options.addConfigOption< bool >("HeatFlux","bool to tell if the wall has constant heat flux (possibly initially), default true.");
  options.addConfigOption< CFuint,Config::DynamicOption<> >("ChangeToIsoT","Iteration after which to switch to an isothermal BC.");
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxMF2D::computeGhostStates(const vector< State* >& intStates,
                                                  vector< State* >& ghostStates,
                                                  const std::vector< RealVector >& normals,
                                                  const std::vector< RealVector >& coords)
{
  const CFuint nbSpecies = m_varSet->getModel()->getNbScalarVars(0);

  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());
  
  ///Maxwell Equations Perfectly Conducting Wall Condition
  const CFreal muZero = 1.2566370614359e-6;
  const CFreal omega = 1.82e5; //14857.14286;
  const CFreal Ez0 = _Ez0;  

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  if (iter >= m_changeToIsoT && m_heatFlux)
  {
    m_heatFlux = false;
  }

  // get some physical data from the model
  const CFreal gamma = m_varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intState   = (*intStates  [iState]);
    State& ghostState = (*ghostStates[iState]);

    //   const CFreal dy = ghostState->getCoordinates()[YY] - innerState->getCoordinates()[YY];
    //  const CFreal dr = MathFunctions::getDistance(ghostState.getCoordinates(),
    //						  intState.getCoordinates());
    //cout<<"dr = "<< dr <<"\n";

    // set the physical data starting from the inner state
    //m_varSet->computePhysicalData(intState,m_intSolPhysData);

    // normal
    const RealVector& normal = normals[iState];

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

  (ghostState)[0] = (intState)[0]; //- dr*muZero*omega*Ez0;	//Bx
  (ghostState)[1] = -(intState)[1]; //By
  (ghostState)[2] = (intState)[2]; //Bz
  (ghostState)[3] = -(intState)[3]; //Ex
  (ghostState)[4] = -(intState)[4]; //Ey
  (ghostState)[5] = -(intState)[5]; //Ez
  (ghostState)[6] = -(intState)[6];	//Psi
  (ghostState)[7] = -(intState)[7];	//Phi

     ///MultiFluidMHD No Slip Isothermal Condition in 2D
     // here a fix is needed in order to have always ghostT > 0
     // if ghostT < 0  then the inner value is set
     const CFuint endEM = 8;
     
    if (m_heatFlux)
    {
 
     //set the densities
     for (CFuint i = 0 ; i < nbSpecies; i++){
       (ghostState)[endEM + i] = (intState)[endEM + i];
     }
 
     //set the Velocities
     for (CFuint i = 0 ; i < nbSpecies; i++){
       (ghostState)[endEM + nbSpecies + 2*i] = -(intState)[endEM + nbSpecies + 2*i];
       (ghostState)[endEM + nbSpecies + 2*i + 1] = -(intState)[endEM + nbSpecies + 2*i + 1];   
     } 
 
     //set the Temperatures
     for (CFuint i = 0 ; i < nbSpecies; i++){
       const CFreal Twall = _wallTemp[i];
       (ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = (intState)[endEM + nbSpecies + 2*nbSpecies + i];
       //cf_assert(2.*Twall - (intState)[endEM + nbSpecies + 2*nbSpecies + i] > 0.); 
     } 
  }
    else
    {
     //set the densities
     for (CFuint i = 0 ; i < nbSpecies; i++){
       (ghostState)[endEM + i] = (intState)[endEM + i];
     }

     //set the Velocities
     for (CFuint i = 0 ; i < nbSpecies; i++){
       (ghostState)[endEM + nbSpecies + 2*i] = -(intState)[endEM + nbSpecies + 2*i];
       (ghostState)[endEM + nbSpecies + 2*i + 1] = -(intState)[endEM + nbSpecies + 2*i + 1]; 
     } 

     //set the Temperatures
     for (CFuint i = 0 ; i < nbSpecies; i++){
       const CFreal Twall = _wallTemp[i];
       //(ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = (intState)[endEM + nbSpecies + 2*nbSpecies + i];
       (ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = 2.*Twall - (intState)[endEM + nbSpecies + 2*nbSpecies + i];
       cf_assert(2.*Twall - (intState)[endEM + nbSpecies + 2*nbSpecies + i] > 0.); 
     } 
    }

    // set the ghost state from its physical data
    //m_varSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxMF2D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
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

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
//    cf_assert(intGrads[iState].size() == 4);

    // normal
    const RealVector& normal = normals[iState];

    // pressure
//     RealVector& pGradI = *intGrads  [iState][0];
//     RealVector& pGradG = *ghostGrads[iState][0];
//     const CFreal nPGrad = pGradI[XX]*normal[XX] + pGradI[YY]*normal[YY];
//     pGradG = pGradI - 2.0*nPGrad*normal;
    *ghostGrads[iState][0] = *intGrads[iState][0];

    // velocity
    *ghostGrads[iState][1] = *intGrads[iState][1];
    *ghostGrads[iState][2] = *intGrads[iState][2];

    if (m_heatFlux)
    {
      // temperature
      RealVector& tempGradI = *intGrads  [iState][3];
      RealVector& tempGradG = *ghostGrads[iState][3];
      const CFreal nTempGrad = tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY];
      //tempGradG = tempGradI - 2.0*nTempGrad * normal +  normal*m_wallQ;
    }
    else
    {
      *ghostGrads[iState][3] = *intGrads[iState][3];
    }
  }
*/
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxMF2D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // get Euler 2D varset
  m_varSet = getMethodData().getUpdateVar().d_castTo< MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  const CFuint nbSpecies = m_varSet->getModel()->getNbScalarVars(0);  

CFLog(VERBOSE, "\n\n\n\n\n\n_wallTemp size = " << _wallTemp.size() << "\n\n\n\n\n");

  for (CFuint i = 0; i < nbSpecies; ++i) {
    cf_assert(_wallTemp[i] >= 0.0);
  }

  if (m_varSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MltiFluid2DMHDVarSet in BCNoSlipWallHeatFluxMF2D!");
  }

  // resize the physical data for internal and ghost solution points
  m_varSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_varSet->getModel()->resizePhysicalData(m_intSolPhysData  );
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

