#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

#include "FluxReconstructionTurb/FluxReconstructionKOmega.hh"
#include "FluxReconstructionTurb/KOmega2DSourceTerm.hh"
#include "KOmega/NavierStokesKOmegaVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::KOmega;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<KOmega2DSourceTerm, FluxReconstructionSolverData, FluxReconstructionKOmegaModule>
KOmega2DSourceTermProvider("KOmega2DSourceTerm");

//////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

KOmega2DSourceTerm::KOmega2DSourceTerm(const std::string& name) :
    StdSourceTerm(name),
    socket_gradients("gradients"),
    socket_wallDistance("wallDistance"),
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
}

//////////////////////////////////////////////////////////////////////////////

KOmega2DSourceTerm::~KOmega2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::computeProductionTerm(const CFuint iState, 
						const CFreal& CoFactor, 
						const CFreal& MUT,
						CFreal& KProdTerm,  
						CFreal& OmegaProdTerm)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Physics::NavierStokes;
  
  SafePtr< NavierStokes2DKOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes2DKOmega >();
  
  const CFuint uID = 1;//getStateVelocityIDs()[XX];
  const CFuint vID = 2;//getStateVelocityIDs()[YY];
  const CFreal dux = (*(m_cellGrads[iState][uID]))[XX];
  const CFreal duy = (*(m_cellGrads[iState][uID]))[YY]; 
  const CFreal dvx = (*(m_cellGrads[iState][vID]))[XX]; 
  const CFreal dvy = (*(m_cellGrads[iState][vID]))[YY]; 
  
  const CFuint nbScalarEqsSets = m_eulerVarSet->getModel()->getNbScalarVarSets();
  const CFuint iK = m_eulerVarSet->getModel()->getFirstScalarVar(nbScalarEqsSets-1);
  const CFreal rho = navierStokesVarSet->getDensity(*((*m_cellStates)[iState]));
  const CFreal avK = m_solPhysData[iK];
  const CFreal coeffTauMu = navierStokesVarSet->getModel().getCoeffTau();
  const CFreal twoThirdRhoK = (2./3.)*(avK * rho);
  
  KProdTerm = coeffTauMu*(MUT*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)+(duy+dvx)*(duy+dvx))))
                             -twoThirdRhoK*(dux+dvy);
  
  ///Production term: Omega
  const CFreal avOmega = m_solPhysData[iK+1];
  const CFreal blendingCoefF1 = navierStokesVarSet->getBlendingCoefficientF1();
  const CFreal sigmaOmega2 = navierStokesVarSet->getSigmaOmega2();
  OmegaProdTerm  = (navierStokesVarSet->getGammaCoef()*rho/MUT) * KProdTerm;
  
  const CFuint kID = (*((*m_cellStates)[iState])).size() - 2;
  const CFuint omegaID = kID + 1;
  
  ///This is used in (BSL,SST), not for normal kOmega
  const CFreal overOmega = 1./avOmega;
  OmegaProdTerm += (1. - blendingCoefF1) * 2. * rho * overOmega * sigmaOmega2*
    MathFunctions::innerProd(*(m_cellGrads[iState][kID]), *(m_cellGrads[iState][omegaID]));
//  OmegaProdTerm *= _Radius; 
  KProdTerm *=CoFactor;
  
  //Make sure negative values dont propagate...
  KProdTerm            = std::max(0., KProdTerm);
  OmegaProdTerm        = std::max(0., OmegaProdTerm);
}
      
////////////////////////////////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::computeDestructionTerm(const CFuint iState, 
						const CFreal& DcoFactor,
						CFreal& K_desterm, 
						CFreal& Omega_desterm)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Physics::NavierStokes;
  
  SafePtr< NavierStokes2DKOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes2DKOmega >();
  
  // check the average state
  const CFuint nbScalarEqsSets = m_eulerVarSet->getModel()->getNbScalarVarSets();
  const CFuint iK = m_eulerVarSet->getModel()->getFirstScalarVar(nbScalarEqsSets-1);
  const CFreal avK     = m_solPhysData[iK];
  const CFreal avOmega = m_solPhysData[iK+1];
  const CFreal rho = navierStokesVarSet->getDensity(*((*m_cellStates)[iState]));
  
  // Destruction term: k
  K_desterm = (-1.) * rho * avOmega * avK * navierStokesVarSet->getBetaStar(*((*m_cellStates)[iState]));
  K_desterm *= DcoFactor; 
  
  // Destruction term: Omega
  Omega_desterm = (-1.) * rho * avOmega * avOmega * navierStokesVarSet->getBeta(*((*m_cellStates)[iState]));
  
  // Make sure negative values dont propagate...
  K_desterm     = std::min(0., K_desterm );
  Omega_desterm = std::min(0., Omega_desterm);
}
      
////////////////////////////////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::addSourceTerm(RealVector& resUpdates)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Physics::NavierStokes; 
  
  const CFuint kID = (*m_cellStates)[0]->size() - 2;
  const CFuint omegaID = kID + 1;
  
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
    
  SafePtr< NavierStokes2DKOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes2DKOmega >();
   
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

    m_currWallDist[iState] = wallDist[stateID];
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
    computeProductionTerm(iSol, 1., mut, m_prodTerm_k, m_prodTerm_Omega);
    computeDestructionTerm(iSol, 1., m_destructionTerm_k, m_destructionTerm_Omega);
    
    m_prodTerm_k     = std::min(10.*fabs(m_destructionTerm_k), m_prodTerm_k);
    //m_prodTerm_Omega = std::min(10.*fabs(m_destructionTerm_Omega), m_prodTerm_Omega);
      
    /// Compute the rhs contribution
    // and Store the unperturbed source terms
    resUpdates[m_nbrEqs*iSol + kID] = m_prodTerm_k + m_destructionTerm_k;
    resUpdates[m_nbrEqs*iSol + omegaID] = m_prodTerm_Omega + m_destructionTerm_Omega;
  }
}

//////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  StdSourceTerm::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::setup()
{
  CFAUTOTRACE;
  StdSourceTerm::setup();

  // get dimensionality
  m_dim = PhysicalModelStack::getActive()->getDim ();
  
  // get NS 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in KOmega2DSourceTerm!");
  }
  cf_assert(m_nbrEqs == 6);
  
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
}

//////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::unsetup()
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
    KOmega2DSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  result.push_back(&socket_gradients);
  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
