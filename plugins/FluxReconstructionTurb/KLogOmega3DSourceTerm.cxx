#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"

#include "FluxReconstructionTurb/FluxReconstructionKOmega.hh"
#include "FluxReconstructionTurb/KLogOmega3DSourceTerm.hh"
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

MethodCommandProvider<KLogOmega3DSourceTerm, FluxReconstructionSolverData, FluxReconstructionKOmegaModule>
KLogOmega3DSourceTermProvider("KLogOmega3DSourceTerm");

//////////////////////////////////////////////////////////////////////////////

void KLogOmega3DSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("LimitProductionTerm","Limit the production terms for stability (Default = True)");
  options.addConfigOption< bool >("BlockDecoupledJacob","Block decouple ST Jacob in NS-KOmega-GammaRe blocks.");
  options.addConfigOption< bool >("IsAxiSym","If axisymmetric, put to true.");
  options.addConfigOption< bool >("IsSSTV","If using SST-V instead of SST-2003, put to true.");
  options.addConfigOption< bool >("NeglectSSTVTerm","If using SST-Vm instead of SST-V, put to true.");
}

//////////////////////////////////////////////////////////////////////////////

KLogOmega3DSourceTerm::KLogOmega3DSourceTerm(const std::string& name) :
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
    m_currWallDist()
{
  addConfigOptionsTo(this);
  
  m_limitP = true;
  setParameter("LimitProductionTerm",&m_limitP);
  
  m_blockDecoupled = false;
  setParameter("BlockDecoupledJacob",&m_blockDecoupled);
  
  m_isAxisymmetric = false;
  setParameter("IsAxiSym",&m_isAxisymmetric);
  
  m_isSSTV = false;
  setParameter("IsSSTV",&m_isSSTV);
  
  m_neglectSSTVTerm = false;
  setParameter("NeglectSSTVTerm",&m_neglectSSTVTerm);
  
  m_overRadius = 1.0;
  
  m_vOverRadius = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

KLogOmega3DSourceTerm::~KLogOmega3DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void KLogOmega3DSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void  KLogOmega3DSourceTerm::getSToStateJacobian(const CFuint iState)
{
  // reset the jacobian
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_stateJacobian[iEq] = 0.0;
  }
  
  SafePtr< NavierStokes3DKLogOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes3DKLogOmega >();
  
  // Set the wall distance before computing the turbulent viscosity
  navierStokesVarSet->setWallDistance(m_currWallDist[iState]);
  
  if (m_isAxisymmetric)
  {
    m_overRadius = 1.0/max(((*m_cellStates)[iState]->getCoordinates())[YY],1.0e-4);
    m_vOverRadius = m_overRadius*(*((*m_cellStates)[iState]))[2];
  }
    
  /// destruction term of k
    
  const CFreal betaStar = navierStokesVarSet->getBetaStar(*((*m_cellStates)[iState]));
  const CFreal avK     = max((*((*m_cellStates)[iState]))[4],0.0);
  const CFreal avOmega = exp((*((*m_cellStates)[iState]))[5]);
  const CFreal T = max(0.0,(*((*m_cellStates)[iState]))[3]);
  const CFreal p = max(0.0,(*((*m_cellStates)[iState]))[0]);
  const CFreal rho = max(navierStokesVarSet->getDensity(*((*m_cellStates)[iState])),0.0);
  const CFreal R = m_eulerVarSet->getModel()->getR();
  const CFreal overRT = 1.0/(R*T);
  const CFreal pOverRTT = p/(R*T*T);
  const CFreal overOmega = 1.0/avOmega;
  const CFreal DkTerm = avOmega * betaStar;
  const CFreal mut = navierStokesVarSet->getTurbDynViscosityFromGradientVars(*((*m_cellStates)[iState]), m_cellGrads[iState]);
  
  const CFreal Dk = -rho * avK * DkTerm;
  
  if (Dk < 0.0)
  {
 
  CFreal tempSTTerm = 0.0;
  
  if (!m_blockDecoupled)
  {
  //p
  tempSTTerm = -avK * overRT * DkTerm;
  m_stateJacobian[0][4] += tempSTTerm;
  m_stateJacobian[0][3] -= tempSTTerm;
  
  //T
  tempSTTerm = DkTerm * avK * pOverRTT;
  m_stateJacobian[3][4] += tempSTTerm;
  m_stateJacobian[3][3] -= tempSTTerm;
  }
  
  //k
  tempSTTerm = -DkTerm * rho;
  m_stateJacobian[4][4] += tempSTTerm;
  m_stateJacobian[4][3] -= tempSTTerm;
  
  //logOmega
  tempSTTerm = -DkTerm * avK * rho;
  m_stateJacobian[5][4] += tempSTTerm;
  m_stateJacobian[5][3] -= tempSTTerm;
  }
  
  /// destruction term of logOmega
  
  const CFreal beta = navierStokesVarSet->getBeta(*((*m_cellStates)[iState]));
  
  const CFreal DomegaTerm = avOmega * beta;
  
  const CFreal Domega = -DomegaTerm * rho;
  
  if (Domega < 0.0)
  {
  if (!m_blockDecoupled)
  {
  //p
  m_stateJacobian[0][5] = -DomegaTerm * overRT;
  
  //T
  m_stateJacobian[3][5] = DomegaTerm * pOverRTT;
  }
  
  //logOmega
  m_stateJacobian[5][5] = -DomegaTerm * rho;
  }
  
  /// production term of k
  
  const CFuint uID = 1;//getStateVelocityIDs()[XX];
  const CFuint vID = 2;//getStateVelocityIDs()[YY];
  const CFreal dux = (*(m_cellGrads[iState][uID]))[XX];
  const CFreal duy = (*(m_cellGrads[iState][uID]))[YY]; 
  const CFreal dvx = (*(m_cellGrads[iState][vID]))[XX]; 
  const CFreal dvy = (*(m_cellGrads[iState][vID]))[YY]; 
  
  const CFreal u = (*((*m_cellStates)[iState]))[1];
  const CFreal v = (*((*m_cellStates)[iState]))[2];
  
  const CFreal gammaIsentropic = m_eulerVarSet->getModel()->getGamma();

  const CFreal coeffTauMu = navierStokesVarSet->getModel().getCoeffTau();
    
  CFreal mutTerm;
  CFreal Pk;
  
  if (!m_isSSTV)
  {
    mutTerm = coeffTauMu*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)-(dux+dvy-m_vOverRadius)*m_vOverRadius)+(duy+dvx)*(duy+dvx));

    Pk = (mutTerm*mut - (2./3.)*(avK * rho)*(dux+dvy+m_vOverRadius));
  }
  else if (m_neglectSSTVTerm)
  {
    mutTerm = coeffTauMu*(duy-dvx)*(duy-dvx);

    Pk = mutTerm*mut;  
  }
  else
  {
    mutTerm = coeffTauMu*(duy-dvx)*(duy-dvx);

    Pk = (mutTerm*mut - (2./3.)*(avK * rho)*(dux+dvy+m_vOverRadius)); 
  }
  
  const CFreal twoThirdduxduy = (2./3.)*(dux+dvy+m_vOverRadius);
  
  if (Pk > 0.0)
  {
  CFreal tempSTTerm = 0.0;
      
  if (!m_blockDecoupled)
  {
  //p
  tempSTTerm = mutTerm*avK*overRT*overOmega;
  if (!m_neglectSSTVTerm) tempSTTerm += -twoThirdduxduy*avK*overRT;
  
  m_stateJacobian[0][4] += tempSTTerm;
  m_stateJacobian[0][3] -= tempSTTerm;
  
  //T
  tempSTTerm = -mutTerm*avK*overOmega*pOverRTT;
  if (!m_neglectSSTVTerm) tempSTTerm += twoThirdduxduy*avK*pOverRTT;
  
  m_stateJacobian[3][4] += tempSTTerm;
  m_stateJacobian[3][3] -= tempSTTerm;
  
  if (m_isAxisymmetric && !m_isSSTV)
      {
        //v
        tempSTTerm = mut*coeffTauMu*4./3.*m_overRadius*(-(dux+dvy)+v*2.0)- (2./3.)*(avK * rho)*m_overRadius;
        m_stateJacobian[2][4] += tempSTTerm;
        m_stateJacobian[2][3] -= tempSTTerm;
      }
  }
  
  //k
  tempSTTerm = mutTerm*rho*overOmega;
  if (!m_neglectSSTVTerm) tempSTTerm += -twoThirdduxduy*rho;
  
  m_stateJacobian[4][4] += tempSTTerm;
  m_stateJacobian[4][3] -= tempSTTerm;
  
  //logOmega
  tempSTTerm = -mutTerm*rho*avK*overOmega;
  
  m_stateJacobian[5][4] += tempSTTerm;
  m_stateJacobian[5][3] -= tempSTTerm;
  }
  
  /// production of logOmega
  
  const CFreal gamma = navierStokesVarSet->getGammaCoef();
  const CFreal blendingCoefF1 = navierStokesVarSet->getBlendingCoefficientF1();
  const CFreal sigmaOmega2 = navierStokesVarSet->getSigmaOmega2();
  
  const CFreal pOmegaFactor = (1. - blendingCoefF1) * 2. * sigmaOmega2* ((*(m_cellGrads[iState][4]))[XX]*(*(m_cellGrads[iState][5]))[XX] + (*(m_cellGrads[iState][4]))[YY]*(*(m_cellGrads[iState][5]))[YY]);
  
  const CFreal Pomega = (gamma*rho/mut) * Pk * overOmega + rho * overOmega * pOmegaFactor;
  
  if (Pomega > 0.0)
  {
  if (!m_blockDecoupled)
  {
  //p
  m_stateJacobian[0][5] += gamma*(mutTerm*overRT*overOmega) + pOmegaFactor*overRT*overOmega;
  if (!m_neglectSSTVTerm) m_stateJacobian[0][5] += -gamma*twoThirdduxduy*overRT;
  
  //T
  m_stateJacobian[3][5] += gamma*pOverRTT*(-mutTerm*overOmega) - pOmegaFactor*pOverRTT*overOmega;
  if (!m_neglectSSTVTerm) m_stateJacobian[3][5] += gamma*pOverRTT*twoThirdduxduy;
  }
 
  //logOmega
  m_stateJacobian[5][5] +=  -gamma*rho*mutTerm*overOmega - pOmegaFactor*rho*overOmega;
  }
  
  if (m_isAxisymmetric)
  { 
      const CFreal avV = sqrt((*((*m_cellStates)[iState]))[1]*(*((*m_cellStates)[iState]))[1]+(*((*m_cellStates)[iState]))[2]*(*((*m_cellStates)[iState]))[2]);

      const CFreal mu = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iState]));
      
      const CFreal rhovr = rho*v*m_overRadius;
      const CFreal vrOverRT = v*m_overRadius*overRT;
      const CFreal vrPOverRTT = v*m_overRadius*pOverRTT;
      const CFreal rhor = m_overRadius*rho;
      
      const CFreal coeffTauMuAS = coeffTauMu*(mu+mut);
      
      //rho
      m_stateJacobian[0][0] += -vrOverRT;
      m_stateJacobian[2][0] += -rhor;
      m_stateJacobian[3][0] += vrPOverRTT;
      
      //rho u
      m_stateJacobian[0][1] += -vrOverRT*u;
      m_stateJacobian[2][1] += -rhor*u;
      m_stateJacobian[3][1] += vrPOverRTT*u;
      m_stateJacobian[1][1] += -rhovr;
      
      //rho v
      m_stateJacobian[0][2] += -vrOverRT*v;
      m_stateJacobian[2][2] += -2.0*rhovr-m_overRadius*m_overRadius*2.0*coeffTauMuAS;
      m_stateJacobian[3][2] += vrPOverRTT*v;
      
      //rho e
      const CFreal isentropicTerm = gammaIsentropic/(gammaIsentropic-1.0);
      m_stateJacobian[0][3] += -isentropicTerm*m_vOverRadius-0.5*avV*avV*vrOverRT;
      m_stateJacobian[2][3] += -isentropicTerm*p*m_overRadius-0.5*(avV*avV+2.0*v*v)*rho*m_overRadius - 4.0/3.0*m_overRadius*m_overRadius*coeffTauMuAS*v;
      m_stateJacobian[3][3] += 0.5*avV*avV*vrPOverRTT;
      m_stateJacobian[1][3] += -rhovr*u;
      
      //rho k
      if (!m_blockDecoupled)
      {
        m_stateJacobian[0][4] += -vrOverRT*avK;
        m_stateJacobian[2][4] += -rhor*avK;
        m_stateJacobian[3][4] += vrPOverRTT*avK;
      }
      m_stateJacobian[4][4] += -rhovr;
      
      //rho logOmega
      const CFreal logOmega = (*((*m_cellStates)[iState]))[5];
      if (!m_blockDecoupled)
      {
        m_stateJacobian[0][5] += -vrOverRT*logOmega;
        m_stateJacobian[2][5] += -rhor*logOmega;
        m_stateJacobian[3][5] += vrPOverRTT*logOmega;
      }
      m_stateJacobian[5][5] += -rhovr;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void  KLogOmega3DSourceTerm::getSToGradJacobian(const CFuint iState)
{
  // reset the jacobian
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_stateJacobian[iEq] = 0.0;
  }
  
  SafePtr< NavierStokes3DKLogOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes3DKLogOmega >();
        
  const CFreal betaStar = navierStokesVarSet->getBetaStar(*((*m_cellStates)[iState]));
  const CFreal avK     = std::max((*((*m_cellStates)[iState]))[4],0.0);
  const CFreal avOmega = std::exp((*((*m_cellStates)[iState]))[5]);
  const CFreal T = max(0.0,(*((*m_cellStates)[iState]))[3]);
  const CFreal p = max(0.0,(*((*m_cellStates)[iState]))[0]);
  const CFreal rho = max(navierStokesVarSet->getDensity(*((*m_cellStates)[iState])),0.0);
  const CFreal R = m_eulerVarSet->getModel()->getR();
  const CFreal overRT = 1.0/(R*T);
  
  /// production term of k
  
  const CFuint uID = 1;//getStateVelocityIDs()[XX];
  const CFuint vID = 2;//getStateVelocityIDs()[YY];
  const CFreal dux = (*(m_cellGrads[iState][uID]))[XX];
  const CFreal duy = (*(m_cellGrads[iState][uID]))[YY]; 
  const CFreal dvx = (*(m_cellGrads[iState][vID]))[XX]; 
  const CFreal dvy = (*(m_cellGrads[iState][vID]))[YY]; 

  const CFreal coeffTauMu = navierStokesVarSet->getModel().getCoeffTau();
  const CFreal twoThirdRhoK = (2./3.)*(avK * rho);
  const CFreal mutTerm = coeffTauMu*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy))+(duy+dvx)*(duy+dvx));
  
  if (!m_blockDecoupled)
  {
  //p
  m_stateJacobian[0][4] += mutTerm*avK*overRT/avOmega - (2./3.)*avK*(dux+dvy)*coeffTauMu*overRT;
  
  //T
  m_stateJacobian[3][4] += -mutTerm*avK*p/(R*T*T*avOmega) + (2./3.)*avK*(dux+dvy)*coeffTauMu*p/(R*T*T);
  }
  
  //k
  m_stateJacobian[4][4] += -(2./3.)*rho*(dux+dvy)*coeffTauMu + mutTerm*rho/avOmega;
  
  //logOmega
  m_stateJacobian[5][4] +=  -mutTerm*rho*avK/avOmega;
  
  /// production of logOmega
  
  const CFreal blendingCoefF1 = navierStokesVarSet->getBlendingCoefficientF1();
  const CFreal sigmaOmega2 = navierStokesVarSet->getSigmaOmega2();
  
  const CFreal pOmegaFactor = (1. - blendingCoefF1) * 2. * sigmaOmega2* ((*(m_cellGrads[iState][4]))[XX]*(*(m_cellGrads[iState][5]))[XX] + (*(m_cellGrads[iState][4]))[YY]*(*(m_cellGrads[iState][5]))[YY]);
  
  if (!m_blockDecoupled)
  {
  //p
  m_stateJacobian[0][5] += navierStokesVarSet->getGammaCoef()*(mutTerm*overRT/avOmega - (2./3.)*(dux+dvy)*coeffTauMu*overRT);
  
  //T
  m_stateJacobian[3][5] += navierStokesVarSet->getGammaCoef()*(-mutTerm*p/(R*T*T*avOmega) + (2./3.)*(dux+dvy)*coeffTauMu*p/(R*T*T));
  }
  
  //k
  //m_stateJacobian[4][5] += 0.0;
  
  //logOmega
  m_stateJacobian[5][5] +=  -navierStokesVarSet->getGammaCoef()*rho*mutTerm/avOmega;
  
  
  
}

//////////////////////////////////////////////////////////////////////////////

void KLogOmega3DSourceTerm::computeProductionTerm(const CFuint iState, 
						const CFreal& CoFactor, 
						const CFreal& MUT,
						CFreal& KProdTerm,  
						CFreal& OmegaProdTerm)
{ 
  SafePtr< NavierStokes3DKLogOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes3DKLogOmega >();
  
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
  
  if (!m_isSSTV)
  {
    KProdTerm = coeffTauMu*(MUT*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)-(dux+dvy-m_vOverRadius)*m_vOverRadius)
			       +(duy+dvx)*(duy+dvx)))
                             -twoThirdRhoK*(dux+dvy+m_vOverRadius);
  }
  else if (m_neglectSSTVTerm)
  {
    const CFreal vorticity2 = (duy-dvx)*(duy-dvx);

    KProdTerm = coeffTauMu*MUT*vorticity2;  
  }
  else
  {
    const CFreal vorticity2 = (duy-dvx)*(duy-dvx);

    KProdTerm = coeffTauMu*MUT*vorticity2 - twoThirdRhoK*(dux+dvy+m_vOverRadius);  
  }
 
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

void KLogOmega3DSourceTerm::computeDestructionTerm(const CFuint iState, 
						const CFreal& DcoFactor,
						CFreal& K_desterm, 
						CFreal& Omega_desterm)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Physics::NavierStokes;
  
  SafePtr< NavierStokes3DKLogOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes3DKLogOmega >();
  
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

void KLogOmega3DSourceTerm::addSourceTerm(RealVector& resUpdates)
{ 
  const CFuint kID = (*m_cellStates)[0]->size() - 2;
  const CFuint omegaID = kID + 1;
  
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
    
  SafePtr< NavierStokes3DKLogOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes3DKLogOmega >();
   
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
    
    if (m_isAxisymmetric)
    {
      m_overRadius = 1.0/max(((*m_cellStates)[iSol]->getCoordinates())[YY],1.0e-4);
      m_vOverRadius = m_overRadius*m_solPhysData[EulerTerm::VY];
    }
      
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
      wallShearStressVelocity[(((*m_cellStates)[iSol]))->getLocalID()] = mut;//sqrt(nuTot*fabs(dUdY));//m_prodTerm_Omega;//std::min(m_prodTerm_Omega,200.0);//m_prodTerm_Omega;//
      
      
      
//      if (m_currWallDist[iSol] < 3.0e-3 && m_currWallDist[iSol] > 2.8e-3)
//  { 
//      CFLog(INFO, "Pk: " << m_prodTerm_k << ", Dk: " << m_destructionTerm_k << ", Pw: " << m_prodTerm_Omega << ", Dw: " << m_destructionTerm_Omega << ", d: " << m_currWallDist[iSol] <<"\n");
//  }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void KLogOmega3DSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  StdSourceTerm::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void KLogOmega3DSourceTerm::setup()
{
  CFAUTOTRACE;
  StdSourceTerm::setup();

  // get dimensionality
  m_dim = PhysicalModelStack::getActive()->getDim ();
  
  // get NS 3D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();
  
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSet in KLogOmega3DSourceTerm!");
  }
  cf_assert(m_nbrEqs == 7 || m_nbrEqs == 9);
  
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

void KLogOmega3DSourceTerm::unsetup()
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
    KLogOmega3DSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  result.push_back(&socket_gradients);
  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  KLogOmega3DSourceTerm::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_wallShearStressVelocity);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
