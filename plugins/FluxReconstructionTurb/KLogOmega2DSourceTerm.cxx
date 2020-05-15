#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

#include "FluxReconstructionTurb/FluxReconstructionKOmega.hh"
#include "FluxReconstructionTurb/KLogOmega2DSourceTerm.hh"
#include "KOmega/NavierStokesKLogOmegaVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::KOmega;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<KLogOmega2DSourceTerm, FluxReconstructionSolverData, FluxReconstructionKOmegaModule>
KLogOmega2DSourceTermProvider("KLogOmega2DSourceTerm");

//////////////////////////////////////////////////////////////////////////////

void KLogOmega2DSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("LimitProductionTerm","Limit the production terms for stability (Default = True)");
}

//////////////////////////////////////////////////////////////////////////////

KLogOmega2DSourceTerm::KLogOmega2DSourceTerm(const std::string& name) :
    StdSourceTerm(name),
    socket_gradients("gradients"),
    socket_wallDistance("wallDistance"),
    socket_wallShearStressVelocity("wallShearStressVelocity"),
    m_dim(),
    m_eulerVarSet(CFNULL),
    m_diffVarSet(CFNULL),
    m_solPhysData(),
    m_dummyGradients(),
    m_cellGrads(),
    m_prodTerm_k(),
    m_prodTerm_Omega(),
    m_destructionTerm_Omega(),
    m_destructionTerm_k(),
    m_currWallDist(),
    m_isAxisymmetric()
{
  addConfigOptionsTo(this);
  
  m_limitP = true;
  setParameter("LimitProductionTerm",&m_limitP);
}

//////////////////////////////////////////////////////////////////////////////

KLogOmega2DSourceTerm::~KLogOmega2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void KLogOmega2DSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void KLogOmega2DSourceTerm::computeProductionTerm(const CFuint iState, 
						const CFreal& CoFactor, 
						const CFreal& MUT,
						CFreal& KProdTerm,  
						CFreal& OmegaProdTerm)
{ 
  SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
  const CFuint uID = 1;//getStateVelocityIDs()[XX];
  const CFuint vID = 2;//getStateVelocityIDs()[YY];
  const CFreal dux = (*(m_cellGrads[iState][uID]))[XX];
  const CFreal duy = (*(m_cellGrads[iState][uID]))[YY]; 
  const CFreal dvx = (*(m_cellGrads[iState][vID]))[XX]; 
  const CFreal dvy = (*(m_cellGrads[iState][vID]))[YY]; 
  
  const CFuint nbScalarEqsSets = m_eulerVarSet->getModel()->getNbScalarVarSets();
  const CFuint iK = m_eulerVarSet->getModel()->getFirstScalarVar(nbScalarEqsSets-1);
  const CFreal rho = navierStokesVarSet->getDensity(*((*m_cellStates)[iState]));
  const CFreal avK = max(m_solPhysData[iK],0.0);
  const CFreal coeffTauMu = navierStokesVarSet->getModel().getCoeffTau();
  const CFreal twoThirdRhoK = (2./3.)*(avK * rho);
  
  KProdTerm = coeffTauMu*(MUT*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy))
			       +(duy+dvx)*(duy+dvx)))
                             -twoThirdRhoK*(dux+dvy);
 
//  KProdTerm = coeffTauMu*(MUT*((4./3.)*(dux*dux + dvy*dvy - dux*dvy) + duy*duy + dvx*dvx + 2*duy*dvx))
//                             -twoThirdRhoK*(dux+dvy);
  
  ///Production term: Omega
  const CFreal avOmega = exp(m_solPhysData[iK+1]);
  const CFreal blendingCoefF1 = navierStokesVarSet->getBlendingCoefficientF1();
  const CFreal sigmaOmega2 = navierStokesVarSet->getSigmaOmega2();
  
  const CFreal overOmega = 1./avOmega;
  OmegaProdTerm  = (navierStokesVarSet->getGammaCoef()*rho/MUT) * KProdTerm * overOmega;
  
  const CFuint kID = (*((*m_cellStates)[iState])).size() - 2;
  const CFuint omegaID = kID + 1;
  
  if (m_limitP)
  {
    KProdTerm     = std::min(10.*fabs(m_destructionTerm_k), KProdTerm);
    OmegaProdTerm = std::min(10.*fabs(m_destructionTerm_k)*overOmega*(navierStokesVarSet->getGammaCoef()*rho/MUT), OmegaProdTerm);
  }
  
  ///This is used in (BSL,SST), not for normal kOmeg
  OmegaProdTerm += (1. - blendingCoefF1) * 2. * rho * overOmega * sigmaOmega2* ((*(m_cellGrads[iState][kID]))[XX]*(*(m_cellGrads[iState][omegaID]))[XX] + (*(m_cellGrads[iState][kID]))[YY]*(*(m_cellGrads[iState][omegaID]))[YY]);
  //OmegaProdTerm += (1. - blendingCoefF1) * 2. * rho * overOmega * sigmaOmega2* ((*(m_cellGrads[iState][kID]))[XX]*(*(m_cellGrads[iState][omegaID]))[XX] + (*(m_cellGrads[iState][kID]))[YY]*(*(m_cellGrads[iState][omegaID]))[YY]);
    //MathFunctions::innerProd(*(m_cellGrads[iState][kID]), *(m_cellGrads[iState][omegaID]));
//  OmegaProdTerm *= _Radius; 
  
  const CFreal mu = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iState]));
   
  const CFreal coeffTauMu1 = coeffTauMu*mu;
  const CFreal coeffTauMu3 = coeffTauMu*MUT*navierStokesVarSet->getSigmaOmega();//(blendingCoefF1 * navierStokesVarSet->getSigmaOmega1() + (1.0-blendingCoefF1)*sigmaOmega2);
  
  OmegaProdTerm += (coeffTauMu1 + coeffTauMu3)*((*(m_cellGrads[iState][omegaID]))[XX]*(*(m_cellGrads[iState][omegaID]))[XX] + (*(m_cellGrads[iState][omegaID]))[YY]*(*(m_cellGrads[iState][omegaID]))[YY]);
  
  KProdTerm *= CoFactor;
  
  //Make sure negative values dont propagate...
  KProdTerm            = std::max(0., KProdTerm);
  OmegaProdTerm        = std::max(0., OmegaProdTerm);
}
      
////////////////////////////////////////////////////////////////////////////////////////////////////////

void KLogOmega2DSourceTerm::computeDestructionTerm(const CFuint iState, 
						const CFreal& DcoFactor,
						CFreal& K_desterm, 
						CFreal& Omega_desterm)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Physics::NavierStokes;
  
  SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
  // check the average state
  const CFuint nbScalarEqsSets = m_eulerVarSet->getModel()->getNbScalarVarSets();
  const CFuint iK = m_eulerVarSet->getModel()->getFirstScalarVar(nbScalarEqsSets-1);
  const CFreal avK     = std::max(m_solPhysData[iK],0.0);
  const CFreal avOmega = std::exp(m_solPhysData[iK+1]);
  const CFreal rho = navierStokesVarSet->getDensity(*((*m_cellStates)[iState]));
  
  // Destruction term: k
  K_desterm = (-1.) * rho * avOmega * avK * navierStokesVarSet->getBetaStar(*((*m_cellStates)[iState]));
  K_desterm *= DcoFactor; 
  
  // Destruction term: Omega
  Omega_desterm = (-1.) * rho * avOmega * navierStokesVarSet->getBeta(*((*m_cellStates)[iState]));//(-1.) * rho * avOmega * avOmega * navierStokesVarSet->getBeta(*((*m_cellStates)[iState]));
  
//  if (m_currWallDist[iState] < 3.0e-3 && m_currWallDist[iState] > 2.8e-3 && !m_isPerturbed)
//  { 
//      CFLog(INFO, "WkB: " << Omega_desterm << ", beta: " << navierStokesVarSet->getBetaStar(*((*m_cellStates)[iState])) << ", rho: " << rho << ", fact: " << DcoFactor << "\n");
//  }
  
  // Make sure negative values dont propagate...
  K_desterm     = std::min(0., K_desterm );
  Omega_desterm = std::min(0., Omega_desterm);
}
      
////////////////////////////////////////////////////////////////////////////////////////////////////////

void KLogOmega2DSourceTerm::addSourceTerm(RealVector& resUpdates)
{ 
  const CFuint kID = (*m_cellStates)[0]->size() - 2;
  const CFuint omegaID = kID + 1;
  
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
    
  SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes2DKLogOmega >();
   
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // set gradients
  const CFuint nbrStates = m_cellStates->size();

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const CFuint stateID = (*m_cellStates)[iState]->getLocalID();
    
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      *(m_cellGrads[iState][iEq]) = gradients[stateID][iEq];
    }
    
    // Get the wall distance
    DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();

    m_currWallDist[iState] = wallDist[stateID];//((*m_cellStates)[iState]->getCoordinates())[YY];
  }
  
  const EquationSubSysDescriptor& eqData = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint nbEqs = eqData.getNbEqsSS();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq(); 
  
  for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
  {
    m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]), m_solPhysData);
      
    // Set the wall distance before computing the turbulent viscosity
    navierStokesVarSet->setWallDistance(m_currWallDist[iSol]);
    
    const CFreal mut = navierStokesVarSet->getTurbDynViscosityFromGradientVars(*((*m_cellStates)[iSol]), m_cellGrads[iSol]);
        
    navierStokesVarSet->computeBlendingCoefFromGradientVars(*((*m_cellStates)[iSol]), *(m_cellGrads[iSol][kID]), *(m_cellGrads[iSol][omegaID]));
 
    //Compute Reynolds stress tensor 
    computeDestructionTerm(iSol, 1., m_destructionTerm_k, m_destructionTerm_Omega);
    computeProductionTerm(iSol, 1., mut, m_prodTerm_k, m_prodTerm_Omega); 
      
    /// Compute the rhs contribution
    // and Store the unperturbed source terms
    resUpdates[m_nbrEqs*iSol + kID] = m_prodTerm_k + m_destructionTerm_k;
    resUpdates[m_nbrEqs*iSol + omegaID] = m_prodTerm_Omega + m_destructionTerm_Omega;
    
    resUpdates[m_nbrEqs*iSol + 3] = -m_prodTerm_k - m_destructionTerm_k;
    
    if (!m_isPerturbed)
    {
      // store the shear stress velocity
      DataHandle< CFreal > wallShearStressVelocity = socket_wallShearStressVelocity.getDataHandle();
      
      const CFreal mu = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iSol]));
      
      const CFreal rho = navierStokesVarSet->getDensity(*((*m_cellStates)[iSol]));
    
      const CFreal nuTot = (mu + mut)/rho;
      
      const CFreal dUdY = (*(m_cellGrads[iSol][1]))[YY];
            
      // take the absolute value of dUdY to avoid nan which causes tecplot to be unable to load the file
      wallShearStressVelocity[(((*m_cellStates)[iSol]))->getLocalID()] = sqrt(nuTot*fabs(dUdY));//m_prodTerm_Omega;//std::min(m_prodTerm_Omega,200.0);//m_prodTerm_Omega;//
      
      
      
//      if (m_currWallDist[iSol] < 3.0e-3 && m_currWallDist[iSol] > 2.8e-3)
//  { 
//      CFLog(INFO, "Pk: " << m_prodTerm_k << ", Dk: " << m_destructionTerm_k << ", Pw: " << m_prodTerm_Omega << ", Dw: " << m_destructionTerm_Omega << ", d: " << m_currWallDist[iSol] <<"\n");
//  }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void KLogOmega2DSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  StdSourceTerm::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void KLogOmega2DSourceTerm::setup()
{
  CFAUTOTRACE;
  StdSourceTerm::setup();

  // get dimensionality
  m_dim = PhysicalModelStack::getActive()->getDim ();
  
  // get NS 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in KLogOmega2DSourceTerm!");
  }
  cf_assert(m_nbrEqs == 6 || m_nbrEqs == 8);
  
  m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData);
  
  m_currWallDist.resize(m_nbrSolPnts);
  
  m_diffVarSet = getMethodData().getDiffusiveVar();
  
  // size cell gradients vector
  m_cellGrads.resize(m_nbrSolPnts);
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    m_cellGrads[iState].resize(m_nbrEqs);
    
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_cellGrads[iState][iEq] = new RealVector(m_dim);
    }
  }
  
  DataHandle< CFreal > wallShearStressVelocity = socket_wallShearStressVelocity.getDataHandle();
  
  DataHandle< CFreal > wallDistance = socket_wallDistance.getDataHandle();
  
  // resize socket
  wallShearStressVelocity.resize(wallDistance.size());
}

//////////////////////////////////////////////////////////////////////////////

void KLogOmega2DSourceTerm::unsetup()
{
  CFAUTOTRACE;
  StdSourceTerm::unsetup();
  
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      deletePtr(m_cellGrads[iState][iVar]); 
    }
    
    m_cellGrads[iState].clear();
  }
  
  m_cellGrads.clear();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    KLogOmega2DSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  result.push_back(&socket_gradients);
  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  KLogOmega2DSourceTerm::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_wallShearStressVelocity);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
