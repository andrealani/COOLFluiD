#include "FluxReconstructionTurb/GammaAlpha3DSourceTerm.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "FluxReconstructionTurb/FluxReconstructionKOmega.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluxReconstructionTurb/KLogOmega3DSourceTerm.hh"
#include "GammaAlpha/NavierStokes3DGammaAlphaPuvt.hh"

#include "Framework/MethodCommandProvider.hh"
#include "MathTools/MathConsts.hh"
#include "KOmega/NavierStokesKLogOmegaVarSetTypes.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::KOmega;
using namespace COOLFluiD::Physics::GammaAlpha;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

////////////////////////////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<GammaAlpha3DSourceTerm, FluxReconstructionSolverData, FluxReconstructionKOmegaModule>
GammaAlpha3DSourceTermFRProvider("GammaAlpha3DSourceTerm");

///////////////////////////////////////////////////////////////////////////////////////////////////

void GammaAlpha3DSourceTerm::defineConfigOptions(Config::OptionList& options)
{
 options.addConfigOption< bool >("PGrad","pressure Gradient");
 options.addConfigOption< bool >("Decouple","Decouple y-ReTheta from k and log(omega), simply solving as fully turbulent.");
 options.addConfigOption< bool >("LimPAlpha","Limit P_alpha.");
 options.addConfigOption< bool >("AddUpdateCoeff","Add the ST time step restriction.");
 options.addConfigOption< bool,Config::DynamicOption<> >("AddDGamma","Add destruction terms for gamma and alpha.");
 options.addConfigOption< CFreal >("LimLambda","Limit Lambda pressure term.");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

GammaAlpha3DSourceTerm::GammaAlpha3DSourceTerm(const std::string& name) :
  KLogOmega3DSourceTerm(name),
  socket_MInfLocal("MInfLocal"),
  socket_uInfLocal("uInfLocal"),
  socket_TInfLocal("TInfLocal"),
  socket_TuInfLocal("TuInfLocal"),
  socket_alphaDiff("alphaDiff"),
  socket_updateCoeff("updateCoeff"),
  m_Rethetat(),
  m_Rethetac(),
  m_Flength(),
  m_vorticity(),
  m_strain(),
  m_order()
{ 
  addConfigOptionsTo(this);
  
  m_PGrad = true;
  setParameter("PGrad",&m_PGrad);
  
  m_decouple = false;
  setParameter("Decouple",&m_decouple);
  
  m_limPRe = false;
  setParameter("LimPAlpha",&m_limPRe);
  
  m_addUpdateCoeff = false;
  setParameter("AddUpdateCoeff",&m_addUpdateCoeff);
  
  m_addDGDA = true;
  setParameter("AddDGamma",&m_addDGDA);
  
  m_lambdaLim = 0.04;
  setParameter("LimLambda",&m_lambdaLim);
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////

GammaAlpha3DSourceTerm::~GammaAlpha3DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void GammaAlpha3DSourceTerm::setup()
{
  KLogOmega3DSourceTerm::setup();
  
  DataHandle< CFreal > MInfLocal = socket_MInfLocal.getDataHandle();
  
  DataHandle< CFreal > uInfLocal = socket_uInfLocal.getDataHandle();
  
  DataHandle< CFreal > TInfLocal = socket_TInfLocal.getDataHandle();
  
  DataHandle< CFreal > TuInfLocal = socket_TuInfLocal.getDataHandle();
  
  DataHandle< CFreal > alphaDiff = socket_alphaDiff.getDataHandle();
  
  DataHandle< CFreal > wallDistance = socket_wallDistance.getDataHandle();
  
  // resize socket
  MInfLocal.resize(wallDistance.size());
  uInfLocal.resize(wallDistance.size());
  TInfLocal.resize(wallDistance.size());
  TuInfLocal.resize(wallDistance.size());
  alphaDiff.resize(wallDistance.size());
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_order = static_cast<CFuint>(order);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void GammaAlpha3DSourceTerm::unsetup()
{
  KLogOmega3DSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  GammaAlpha3DSourceTerm::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = KLogOmega3DSourceTerm::providesSockets();
  result.push_back(&socket_MInfLocal);
  result.push_back(&socket_uInfLocal);
  result.push_back(&socket_TInfLocal);
  result.push_back(&socket_TuInfLocal);
  result.push_back(&socket_alphaDiff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
GammaAlpha3DSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = KLogOmega3DSourceTerm::needsSockets();
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void  GammaAlpha3DSourceTerm::getStrain(const CFuint iState)
{
  const CFreal gradU_X = (*(m_cellGrads[iState][1]))[XX];
  const CFreal gradU_Y = (*(m_cellGrads[iState][1]))[YY];
  const CFreal gradU_Z = (*(m_cellGrads[iState][1]))[ZZ];
  const CFreal gradV_X = (*(m_cellGrads[iState][2]))[XX];
  const CFreal gradV_Y = (*(m_cellGrads[iState][2]))[YY];
  const CFreal gradV_Z = (*(m_cellGrads[iState][2]))[ZZ];
  const CFreal gradW_X = (*(m_cellGrads[iState][3]))[XX];
  const CFreal gradW_Y = (*(m_cellGrads[iState][3]))[YY];
  const CFreal gradW_Z = (*(m_cellGrads[iState][3]))[ZZ];
      
  const CFreal gradSumXY = (gradU_Y + gradV_X);
  const CFreal gradSumXZ = (gradU_Z + gradW_X);
  const CFreal gradSumYZ = (gradW_Y + gradV_Z);
     
  m_strain = std::sqrt(2.*(gradU_X*gradU_X + gradV_Y*gradV_Y + gradW_Z*gradW_Z) + gradSumXY*gradSumXY + gradSumXZ*gradSumXZ + gradSumYZ*gradSumYZ);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void  GammaAlpha3DSourceTerm::getVorticity(const CFuint iState)
{
  const CFreal gradU_Y = (*(m_cellGrads[iState][1]))[YY];
  const CFreal gradU_Z = (*(m_cellGrads[iState][1]))[ZZ];
  const CFreal gradV_X = (*(m_cellGrads[iState][2]))[XX];
  const CFreal gradV_Z = (*(m_cellGrads[iState][2]))[ZZ];
  const CFreal gradW_X = (*(m_cellGrads[iState][3]))[XX];
  const CFreal gradW_Y = (*(m_cellGrads[iState][3]))[YY];
  
  const CFreal gradDiffXY = (gradU_Y - gradV_X);
  const CFreal gradDiffXZ = (gradU_Z - gradW_X);
  const CFreal gradDiffYZ = (gradW_Y - gradV_Z);
  
  m_vorticity = sqrt(gradDiffXY*gradDiffXY+gradDiffXZ*gradDiffXZ+gradDiffYZ*gradDiffYZ);
}

//////////////////////////////////////////////////////////////////////////////

void GammaAlpha3DSourceTerm::addSourceTerm(RealVector& resUpdates)
{       
  const CFuint kID = (*m_cellStates)[0]->size() - 4;
  const CFuint omegaID = kID + 1;
  const CFuint gammaID = kID + 2;
  const CFuint alphaID = kID + 3;

  const CFuint iKPD = m_eulerVarSet->getModel()->getFirstScalarVar(0);
      
  SafePtr< NavierStokes3DGammaAlphaPuvt > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes3DGammaAlphaPuvt >();
   
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
    const CFreal mu = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iSol]));
    
    navierStokesVarSet->computeBlendingCoefFromGradientVars(*((*m_cellStates)[iSol]), *(m_cellGrads[iSol][kID]), *(m_cellGrads[iSol][omegaID]));
    
    // Get Vorticity
    getVorticity(iSol);
    
    // get current state values
    const CFreal avP     = std::max(m_solPhysData[EulerTerm::P],1.0e-3);
    const CFreal avV     = m_solPhysData[EulerTerm::V];
    const CFreal avK     = std::max(m_solPhysData[iKPD],0.0);
    const CFreal avOmega = std::exp(m_solPhysData[iKPD+1]);
    
    const CFreal avGa    = std::min(std::max(m_solPhysData[iKPD+2],0.01),0.99);
    //const CFreal avGa    = std::min(std::max(m_solPhysData[iKPD+2],0.0),1.0);

    const CFreal rho = navierStokesVarSet->getDensity(*((*m_cellStates)[iSol]));
    const CFreal u = m_solPhysData[EulerTerm::VX];
    const CFreal v = m_solPhysData[EulerTerm::VY];
    const CFreal w = m_solPhysData[EulerTerm::VZ];
    const CFreal dpx = (*(m_cellGrads[iSol][0]))[XX];
    const CFreal dpy = (*(m_cellGrads[iSol][0]))[YY];
    const CFreal dpz = (*(m_cellGrads[iSol][0]))[ZZ];
    
    // compute local freestream values with isentropic assumption
    const CFreal MInf = m_eulerVarSet->getModel()->getMachInf();
    const CFreal pInf = m_eulerVarSet->getModel()->getStaticPressInf();//getPressInf();
    const CFreal gammaIsentropic = m_eulerVarSet->getModel()->getGamma();
    const CFreal uInf = m_eulerVarSet->getModel()->getVelInf();
    const CFreal rhoInf = (MInf*MInf)/(uInf*uInf)*gammaIsentropic*pInf;
  
    const CFreal pTotTerm = 1.0+(gammaIsentropic-1.0)/2.0*MInf*MInf;
    const CFreal pTotExponent = gammaIsentropic/(gammaIsentropic-1.0);
    const CFreal pTotalInf = pow(pTotTerm,pTotExponent)*avP;
    
    const CFreal R = m_eulerVarSet->getModel()->getR();
  
    const CFreal rhoInfLocal = rhoInf*pow(avP/pInf,1.0/gammaIsentropic);
    const CFreal MInfLocal = sqrt(2.0/(gammaIsentropic-1.0)*(pow(pTotalInf/avP,1.0/pTotExponent)-1.0));
    const CFreal uInfLocal = sqrt(gammaIsentropic*avP/rhoInfLocal)*MInfLocal;
    const CFreal TInfLocal = avP/(rhoInfLocal*R);
    
    const CFreal Tu = min(max(100.0 * (std::sqrt(2.0*avK/3.0))/(avV),0.027),13.0);
    
    if (!(Tu >= 0.0)) CFLog(INFO, "Tu problem, k: " << avK << ", V: " << avV << "\n");
    
    const CFreal TLocal = (*((*m_cellStates)[iSol]))[4];
    (*((*m_cellStates)[iSol]))[4] = TInfLocal;
    
    const CFreal muInfLocal = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iSol]));
    
    (*((*m_cellStates)[iSol]))[4] = TLocal;
    
    const CFreal alphaMin = rhoInfLocal*uInfLocal*uInfLocal*0.2247/(std::sqrt(rhoInfLocal*muInfLocal*(1.0+0.38*pow(MInfLocal,0.6)))*4900.0);
    
    const CFreal avAlpha    = std::max(m_solPhysData[iKPD+3],alphaMin);    
    
    CFreal ReThetat;

    if (!m_PGrad)
    {
      ReThetat = rhoInfLocal*uInfLocal*uInfLocal*0.2247141166/(avAlpha*sqrt(rhoInfLocal*muInfLocal));
    }
    else
    {
      ReThetat = getRethetatwithPressureGradient(avAlpha,rhoInfLocal,uInfLocal,muInfLocal,mu,rho,dpx,dpy,dpz);
    }

    const CFreal TuInfLocal = getTuInfLocal(ReThetat,MInfLocal,Tu);

    ///@todo
    CFreal fOnset = 1.0;
    
    // compute destruction terms for k and logOmega
    // Destruction term: k
    m_destructionTerm_k = (-1.) * rho * avOmega * avK * navierStokesVarSet->getBetaStar(*((*m_cellStates)[iSol]));
      
    // Destruction term: Omega
    m_destructionTerm_Omega = (-1.) * rho * avOmega * navierStokesVarSet->getBeta(*((*m_cellStates)[iSol]));//(-1.) * rho * avOmega * avOmega * navierStokesVarSet->getBeta(*((*m_cellStates)[iState]));

  
    // Make sure negative values dont propagate...
    m_destructionTerm_k     = std::min(0., m_destructionTerm_k );
    m_destructionTerm_Omega = std::min(0., m_destructionTerm_Omega);
        
    ///The Modified Production  terms for k and logOmega
    const CFuint uID = 1;//getStateVelocityIDs()[XX];
    const CFuint vID = 2;//getStateVelocityIDs()[YY];
    const CFreal dux = (*(m_cellGrads[iSol][uID]))[XX];
    const CFreal duy = (*(m_cellGrads[iSol][uID]))[YY];
    const CFreal duz = (*(m_cellGrads[iSol][uID]))[ZZ];
    const CFreal dvx = (*(m_cellGrads[iSol][vID]))[XX]; 
    const CFreal dvy = (*(m_cellGrads[iSol][vID]))[YY]; 
    const CFreal dvz = (*(m_cellGrads[iSol][vID]))[ZZ];
    const CFreal dwx = (*(m_cellGrads[iSol][3]))[XX]; 
    const CFreal dwy = (*(m_cellGrads[iSol][3]))[YY]; 
    const CFreal dwz = (*(m_cellGrads[iSol][3]))[ZZ];
    const CFreal dgammax = (*(m_cellGrads[iSol][7]))[XX];
    const CFreal dgammay = (*(m_cellGrads[iSol][7]))[YY];
    const CFreal dgammaz = (*(m_cellGrads[iSol][7]))[ZZ];
    const CFreal dax = (*(m_cellGrads[iSol][8]))[XX];
    const CFreal day = (*(m_cellGrads[iSol][8]))[YY];
    const CFreal daz = (*(m_cellGrads[iSol][8]))[ZZ];
    const CFreal dkx = (*(m_cellGrads[iSol][5]))[XX];
    const CFreal dky = (*(m_cellGrads[iSol][5]))[YY];
    const CFreal dkz = (*(m_cellGrads[iSol][5]))[ZZ];

    const CFreal coeffTauMu = navierStokesVarSet->getModel().getCoeffTau();
    const CFreal gammaTerm = avGa+avGa*(1.0-avGa);
    const CFreal twoThirdRhoK = (2./3.)*(avK * rho * gammaTerm);
    
//    if (!m_isSSTV)
//    {
//      m_prodTerm_k = coeffTauMu*(mut*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)-(dux+dvy-m_vOverRadius)*m_vOverRadius)
//			         +(duy+dvx)*(duy+dvx)))
//                                 -twoThirdRhoK*(dux+dvy+m_vOverRadius);
//    }
    if (m_neglectSSTVTerm)
    {
      //const CFreal vorticity2 = (duy-dvx)*(duy-dvx);

      m_prodTerm_k = coeffTauMu*mut*m_vorticity*m_vorticity;  
    }
    else
    {
      //const CFreal vorticity2 = (duy-dvx)*(duy-dvx);

      m_prodTerm_k = coeffTauMu*mut*m_vorticity*m_vorticity - twoThirdRhoK*(dux+dvy+dwz);  
    }
  
//    m_prodTerm_k = coeffTauMu*(mut*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)-(dux+dvy-m_vOverRadius)*m_vOverRadius)
//			         +(duy+dvx)*(duy+dvx)))
//                                 -twoThirdRhoK*(dux+dvy+m_vOverRadius);
  
    ///Production term: Omega
    const CFreal blendingCoefF1 = navierStokesVarSet->getBlendingCoefficientF1();
    const CFreal sigmaOmega2 = navierStokesVarSet->getSigmaOmega2();
     
    const CFreal overOmega = 1./avOmega;
    m_prodTerm_Omega  = (navierStokesVarSet->getGammaCoef()*rho/mut) * m_prodTerm_k * overOmega;
  
    if (m_limitP)
    {
      m_prodTerm_k     = std::min(10.*fabs(m_destructionTerm_k), m_prodTerm_k);
      m_prodTerm_Omega = std::min(10.*fabs(m_destructionTerm_k)*overOmega*(navierStokesVarSet->getGammaCoef()*rho/mut), m_prodTerm_Omega);
    }  
  
    ///This is used in (BSL,SST), not for normal kOmeg
    m_prodTerm_Omega += (1. - blendingCoefF1) * 2. * rho * overOmega * sigmaOmega2* ((*(m_cellGrads[iSol][kID]))[XX]*(*(m_cellGrads[iSol][omegaID]))[XX] + (*(m_cellGrads[iSol][kID]))[YY]*(*(m_cellGrads[iSol][omegaID]))[YY] + (*(m_cellGrads[iSol][kID]))[ZZ]*(*(m_cellGrads[iSol][omegaID]))[ZZ]);
  //OmegaProdTerm += (1. - blendingCoefF1) * 2. * rho * overOmega * sigmaOmega2* ((*(m_cellGrads[iState][kID]))[XX]*(*(m_cellGrads[iState][omegaID]))[XX] + (*(m_cellGrads[iState][kID]))[YY]*(*(m_cellGrads[iState][omegaID]))[YY]);
    //MathFunctions::innerProd(*(m_cellGrads[iState][kID]), *(m_cellGrads[iState][omegaID]));
//  OmegaProdTerm *= _Radius; 
     
    const CFreal coeffTauMu1 = coeffTauMu*mu;
    const CFreal coeffTauMu3 = coeffTauMu*mut*navierStokesVarSet->getSigmaOmega();//(blendingCoefF1 * navierStokesVarSet->getSigmaOmega1() + (1.0-blendingCoefF1)*sigmaOmega2);
  
    m_prodTerm_Omega += (coeffTauMu1 + coeffTauMu3)*((*(m_cellGrads[iSol][omegaID]))[XX]*(*(m_cellGrads[iSol][omegaID]))[XX] + (*(m_cellGrads[iSol][omegaID]))[YY]*(*(m_cellGrads[iSol][omegaID]))[YY] + (*(m_cellGrads[iSol][omegaID]))[ZZ]*(*(m_cellGrads[iSol][omegaID]))[ZZ]);
  
    //Make sure negative values dont propagate...
    m_prodTerm_k            = std::max(0., m_prodTerm_k);
    m_prodTerm_Omega        = std::max(0., m_prodTerm_Omega);
    
    // compute production term of gamma
    const CFreal cfg1 = 1.735;
    const CFreal cfg2 = 5.45;
    const CFreal cfg3 = 0.95375;
    const CFreal cfg4 = 2.2;
    
    CFreal fg = 1.0;
    
    if (avGa < 0.45) fg -= exp(-cfg1*tan(cfg2*avGa-cfg3)-cfg4);
    
    const CFreal fMnsigmaTerm = 1.0+0.58*pow(MInfLocal,0.6);
    const CFreal fMnsigma = 1.0/(fMnsigmaTerm*fMnsigmaTerm);
    
    CFreal fk = 1.0;
    
    CFreal Rethetac = getRethetat(Tu, false);
    
    ///@todo check if this is local M or local MInf
    Rethetac *= sqrt(1.0+0.38*pow(MInfLocal,0.6));
    
    CFreal lambda = 0.0;

    CFreal dpds = 1.0/avV*(u*dpx + v*dpy + w*dpz);
    
    if (m_PGrad)
    {
      const CFreal lambdaTerm1 = -mu/(rho*rho*uInfLocal*uInfLocal*uInfLocal);

      const CFreal dpdsLimit = -m_lambdaLim/(lambdaTerm1*Rethetac*Rethetac);
      dpds = min(max(dpds,-dpdsLimit),dpdsLimit);
    
      const CFreal kPGrad = -muInfLocal/(rhoInfLocal*rhoInfLocal*uInfLocal*uInfLocal*uInfLocal)*fabs(1.0-MInfLocal*MInfLocal)*dpds;
    
      CFreal PRC;
      if (kPGrad < 0)
      {
        const CFreal PRCfactor1 = 474.0*pow(TuInfLocal,-2.9);
        const CFreal PRCexponent = 1.0-exp(2.0e6*kPGrad);
        PRC = pow(PRCfactor1,PRCexponent);
      }
      else
      {
        const CFreal PRCexponent = -3227.0*pow(kPGrad,0.5985);
        PRC = pow(10.0,PRCexponent);
      }
    
      fk += fg*(PRC-1.0); 
      
      lambda = -Rethetac*Rethetac*mu/(rho*rho*uInfLocal*uInfLocal*uInfLocal)*dpds;
    }
    
    const CFreal nsigma = 1.25e-11*pow(TuInfLocal,7.0/4.0)*fk*fMnsigma;
    
    const CFreal beta = sqrt(nsigma)*uInfLocal*rhoInfLocal/muInfLocal;
    
    //const CFreal limGa =  std::min(avGa,0.9999);
    CFreal prodTerm_Ga = fOnset*2.0*fg*(1.0-avGa)*sqrt(-log(1.0-avGa))*beta*rho*avV;
    
    // compute dissipation term of gamma
    const CFreal cEg = 20.0;
  
    const CFreal fMuGamma = 1.0-exp(-256.0*(m_currWallDist[iSol]*uInfLocal*rhoInfLocal/muInfLocal)*(m_currWallDist[iSol]*uInfLocal*rhoInfLocal/muInfLocal));
    const CFreal fMMuGamma = (1.0+0.26*(gammaIsentropic-1.0)/2.0*MInfLocal*MInfLocal)*sqrt(1+0.38*pow(MInfLocal,0.6));
  
    //const CFreal gammaLim = std::min(std::max(0.01,avGa),0.99);
    const CFreal muGamma = 0.57*pow(-log(1.0-avGa),-5.0/6.0*(1.0-avGa))*fMuGamma*fMMuGamma*mu;
    
    //const CFreal dudn = -1.0/(avV*avV)*(u*u*duy-v*v*dvx+u*v*(dvy-dux));//-duy; //
    //CFreal dgammadn = 1.0/avV*(v*dgammax - u*dgammay);//-dgammay; //
    
    const CFreal overVel = 1.0/avV;
    const CFreal dVeldx = overVel*(u*dux+v*dvx+w*dwx);
    const CFreal dVeldy = overVel*(u*duy+v*dvy+w*dwy);
    const CFreal dVeldz = overVel*(u*duz+v*dvz+w*dwz);
    
//    const CFreal dVeldx = overVel*(u*dux+v*dvx);
//    const CFreal dVeldy = overVel*(u*duy+v*dvy);
//    const CFreal dVeldz = overVel*(u*duz+v*dvz);
    
    CFreal dVeldnx;
    CFreal dVeldny;
    CFreal dVeldnz;
    CFreal dVeldnSize;
    
    // temporary solution for cone case, code below based on largest du/dn, i.e. direction for shear vector, should work, but u grads in 3D seem too oscillatory in P1
    const CFreal sinA = 0.139173101;//0.1218693434;    
const CFreal cosA = 0.9902680687;//0.9925461516;
const CFreal xvec = -sinA;
    const CFreal zcurr = ((*m_cellStates)[iSol]->getCoordinates())[ZZ];
    const CFreal ycurr = ((*m_cellStates)[iSol]->getCoordinates())[YY];
    const CFreal nsizecurr = 1.0/sqrt(zcurr*zcurr+ycurr*ycurr);
    const CFreal dVeldn1 = (dVeldy*ycurr*cosA+zcurr*dVeldz*cosA)*nsizecurr+xvec*dVeldx;
    dVeldnSize = fabs(dVeldn1);
    dVeldnx = xvec;
    dVeldny = ycurr*nsizecurr*cosA;//dVeldn1*ycurr*nsizecurr/dVeldnSize;
    dVeldnz = zcurr*nsizecurr*cosA;//dVeldn1*zcurr*nsizecurr/dVeldnSize;
        
//    if (fabs(w) <= fabs(u))
//    {
//        const CFreal nsize1 = 1.0/sqrt(u*u+v*v);
//        const CFreal dVeldn1 = (dVeldx*v-u*dVeldy)*nsize1;
//        
//        const CFreal nsize2 = 1.0/sqrt(w*u*w*u+v*w*v*w+(u*u+v*v)*(u*u+v*v));
//        const CFreal dVeldn2 = (dVeldx*w*u+v*w*dVeldy-(u*u+v*v)*dVeldz)*nsize2;
//        
//        dVeldnSize = sqrt(dVeldn1*dVeldn1+dVeldn2*dVeldn2);
//        
//        dVeldnx = (dVeldn1*v*nsize1+dVeldn2*u*w*nsize2)/dVeldnSize;
//        
//        dVeldny = (-dVeldn1*u*nsize1+dVeldn2*v*w*nsize2)/dVeldnSize;
//        
//        dVeldnz = (-dVeldn2*(u*u+v*v)*nsize2)/dVeldnSize;
//        
////      const CFreal dgammadn1 = 1.0/avV*(v*dgammax - u*dgammay);//-dgammay; //
////      
////      const CFreal dgammadn2 = 1.0/(avV*avV)*(u*w*dgammax + v*w*dgammay - (u*u+v*v)*dgammaz);
////      
////      dgammadn = sqrt(dgammadn1*dgammadn1+dgammadn2*dgammadn2);
//    }
//    else
//    {
//        const CFreal nsize1 = 1.0/sqrt(w*w+v*v);
//        const CFreal dVeldn1 = nsize1*(-dVeldy*w+v*dVeldz);
//        
//        const CFreal nsize2 = 1.0/sqrt(v*u*v*u+u*w*u*w+(w*w+v*v)*(w*w+v*v));
//        const CFreal dVeldn2 = nsize2*(-dVeldy*v*u-u*w*dVeldz+(w*w+v*v)*dVeldx);
//        
//        dVeldnSize = sqrt(dVeldn1*dVeldn1+dVeldn2*dVeldn2);
//        
//        dVeldnx = (dVeldn2*(w*w+v*v)*nsize2)/dVeldnSize;
//        
//        dVeldny = (-dVeldn1*w*nsize1-dVeldn2*v*u*nsize2)/dVeldnSize;
//        
//        dVeldnz = (dVeldn1*v*nsize1-dVeldn2*u*w*nsize2)/dVeldnSize;
//        
////      const CFreal dgammadn1 = 1.0/avV*(-w*dgammax - u*dgammay);//-dgammay; //
////      
////      const CFreal dgammadn2 = 1.0/(avV*avV)*(u*w*dgammax + v*w*dgammay - (u*u+v*v)*dgammaz);
////      
////      dgammadn = sqrt(dgammadn1*dgammadn1+dgammadn2*dgammadn2);
////      
////      const CFreal tauT1 = muTot*(-m_unitNormalFlxPnts[iFlxPnt][ZZ]*MathFunctions::innerProd(m_vGrad,m_unitNormalFlxPnts[iFlxPnt]) +
////                           m_unitNormalFlxPnts[iFlxPnt][YY]*MathFunctions::innerProd(m_wGrad,m_unitNormalFlxPnts[iFlxPnt]));
////      
////      const CFreal tauT2 = muTot*(-m_unitNormalFlxPnts[iFlxPnt][XX]*m_unitNormalFlxPnts[iFlxPnt][YY]*MathFunctions::innerProd(m_vGrad,m_unitNormalFlxPnts[iFlxPnt]) -
////                           m_unitNormalFlxPnts[iFlxPnt][XX]*m_unitNormalFlxPnts[iFlxPnt][ZZ]*MathFunctions::innerProd(m_wGrad,m_unitNormalFlxPnts[iFlxPnt]) + 
////                           (m_unitNormalFlxPnts[iFlxPnt][YY]*m_unitNormalFlxPnts[iFlxPnt][YY]+m_unitNormalFlxPnts[iFlxPnt][ZZ]*m_unitNormalFlxPnts[iFlxPnt][ZZ])*MathFunctions::innerProd(m_uGrad,m_unitNormalFlxPnts[iFlxPnt]));
////      tau = sqrt(tauT1*tauT1+tauT2*tauT2);
//    }
    
    const CFreal dudn = -dVeldnSize;//1.0/(avV*avV)*(u*u*duy-v*v*dvx+u*v*(dvy-dux));//-duy; //
    
    const CFreal dgammadn = (dVeldnSize < 1.0e-7*avV) ? 0.0 : -(dVeldnx*dgammax + dVeldny*dgammay + dVeldnz*dgammaz);//-dgammay; //
    
    //CFLog(INFO, "nx: " << dVeldnx << ", ny: " << dVeldny << ", nz: " << dVeldnz << ", dudn: " << dudn << ", dgdn: " << dgammadn << "\n");
    
    CFreal  destructionTerm_Ga  = (-1.0) * cEg/0.57 * muGamma * avV/(uInfLocal*uInfLocal) * dudn * dgammadn;
    
    //destructionTerm_Ga = min(max(-10.0*fabs(prodTerm_Ga),destructionTerm_Ga),10.0*fabs(prodTerm_Ga));
    //if (avGa<0.1) destructionTerm_Ga = min(max(-fabs(prodTerm_Ga),destructionTerm_Ga),fabs(prodTerm_Ga));//if (avGa<0.4) destructionTerm_Ga = max(-fabs(prodTerm_Ga),destructionTerm_Ga);
    //destructionTerm_Ga = min(0.0,destructionTerm_Ga);
    
    destructionTerm_Ga = max(-17.0*fabs(prodTerm_Ga),destructionTerm_Ga);
    
    // compute production term of alpha
    const CFreal cpa1 = 0.03;
    //const CFreal cpa2 = 100.0;//50.0;
    
    const CFreal alphac = rhoInfLocal*uInfLocal*uInfLocal*pow(0.09+lambda,0.62)/(Rethetac*sqrt(rhoInfLocal*muInfLocal));
    
    const CFreal t = (500.0 * mu )/(rho * avV * avV);
    
    //const CFreal dkdn = 1.0/avV*(v*dkx - u*dky);//-dky; // 
    const CFreal dkdn = (dVeldnSize < 1.0e-7*avV) ? 0.0 : -(dVeldnx*dkx + dVeldny*dky + dVeldnz*dkz);
    
    const CFreal  Rew         = (rho * m_currWallDist[iSol] * m_currWallDist[iSol] * avOmega)/(mu);   
    const CFreal  Fwake1      = (1.0e-5 * Rew)*(1.0e-5 * Rew);   
    const CFreal  Fwake       = exp(-Fwake1);
    const CFreal  thetaBL     = (ReThetat*mu)/(rho*avV);
     
    const CFreal  delta       = (375.0 * m_vorticity * m_currWallDist[iSol] * thetaBL)/(avV);
    
    const CFreal  coefFtheta0 = (m_currWallDist[iSol]/delta)*(m_currWallDist[iSol]/delta)*(m_currWallDist[iSol]/delta)*(m_currWallDist[iSol]/delta);
    const CFreal  coefFtheta1 = exp(-coefFtheta0);
    const CFreal  Ftheta1     = Fwake * coefFtheta1;
    const CFreal  Ftheta3     = 1.0-(((avGa-0.01)/(0.99-0.01))*((avGa-0.01)/(0.99-0.01)));
    const CFreal  Ftheta4     = std::max(Ftheta1,Ftheta3);
    CFreal  Fthetat     = std::min(Ftheta4,1.0);
    
    const CFreal dkdnFactor = fabs(dkdn)/max(avK,1.0e-6);

    if (dkdnFactor>100.0) Fthetat = 1.0;
    
    CFreal prodTerm_alpha = cpa1 * (rho/t) * (alphac - avAlpha) * (1.0 - Fthetat);

    // compute destruction term of alpha
    const CFreal cna = 0.4;
    //const CFreal dadn = 1.0/avV*(v*dax - u*day);
    const CFreal dadn = (dVeldnSize < 1.0e-7*avV) ? 0.0 : -(dVeldnx*dax + dVeldny*day + dVeldnz*daz);
    
    CFreal destructionTerm_alpha = -1.0*cna*rho*avV*dadn;
    
    //destructionTerm_alpha = max(-10.0*fabs(prodTerm_alpha),destructionTerm_alpha);

    ///Make sure negative values dont propagate
    prodTerm_Ga        = max(0., prodTerm_Ga);

    if (m_limPRe) prodTerm_alpha     = max(0., prodTerm_alpha);
    
    if (m_limPRe) destructionTerm_Ga = min(0., destructionTerm_Ga);
    
    if (m_limPRe) destructionTerm_alpha = min(0., destructionTerm_alpha);
      
    /// Compute the rhs contribution
    // and Store the unperturbed source terms
    resUpdates[m_nbrEqs*iSol + kID] = m_prodTerm_k + m_destructionTerm_k;
    resUpdates[m_nbrEqs*iSol + omegaID] = m_prodTerm_Omega + m_destructionTerm_Omega;
    resUpdates[m_nbrEqs*iSol + gammaID] = prodTerm_Ga;// + destructionTerm_Ga;
    resUpdates[m_nbrEqs*iSol + alphaID] = prodTerm_alpha + destructionTerm_alpha;
    
    if (m_addDGDA)
    {
      resUpdates[m_nbrEqs*iSol + gammaID] += destructionTerm_Ga;
    }
    
    resUpdates[m_nbrEqs*iSol + 4] = -m_prodTerm_k - m_destructionTerm_k;
    
    //if (m_solPhysData[iKPD+2]<0) CFLog(INFO, "ST update: " << resUpdates << "\n");
    
    if (!m_isPerturbed)
    {
      // store the shear stress velocity
      DataHandle< CFreal > wallShearStressVelocity = socket_wallShearStressVelocity.getDataHandle();
      
      DataHandle< CFreal > MInfLocalSocket = socket_MInfLocal.getDataHandle();
      
      DataHandle< CFreal > uInfLocalSocket = socket_uInfLocal.getDataHandle();
      
      DataHandle< CFreal > TInfLocalSocket = socket_TInfLocal.getDataHandle();
      
      DataHandle< CFreal > TuInfLocalSocket = socket_TuInfLocal.getDataHandle();
      
      DataHandle< CFreal > alphaDiffSocket = socket_alphaDiff.getDataHandle();
      
      const CFreal mu = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iSol]));
      
      const CFreal rho = navierStokesVarSet->getDensity(*((*m_cellStates)[iSol]));
    
      const CFreal nuTot = (mu + mut)/rho;
      
      const CFreal dUdY = (*(m_cellGrads[iSol][1]))[YY];
            
      // take the absolute value of dUdY to avoid nan which causes tecplot to be unable to load the file
      wallShearStressVelocity[(((*m_cellStates)[iSol]))->getLocalID()] = sqrt(nuTot*fabs(dUdY));//m_prodTerm_Omega;//std::min(m_prodTerm_Omega,200.0);//m_prodTerm_Omega;//
      
      MInfLocalSocket[(((*m_cellStates)[iSol]))->getLocalID()] = dkdnFactor;//MInfLocal;//destructionTerm_Ga;//prodTerm_Ga;//dudn;//
      
      uInfLocalSocket[(((*m_cellStates)[iSol]))->getLocalID()] = 1.0/avV*(u*dpx + v*dpy + w*dpz);//uInfLocal;//muGamma;//destructionTerm_Ga;//dgammadn;//
      
      TInfLocalSocket[(((*m_cellStates)[iSol]))->getLocalID()] = dVeldnx;//dpds;//TInfLocal;//dgammadn;//dudn;//prodTerm_alpha;//dadn;//
      
      TuInfLocalSocket[(((*m_cellStates)[iSol]))->getLocalID()] = dVeldny;//TuInfLocal;//avV/(uInfLocal*uInfLocal);//destructionTerm_alpha;//destructionTerm_Ga;//
      
      alphaDiffSocket[(((*m_cellStates)[iSol]))->getLocalID()] = dVeldnz;//((*m_cellStates)[iSol]->getCoordinates())[YY];//m_currWallDist[iSol];//alphac - avAlpha;//dudn;//muGamma;//
    }
  }
}
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

CFreal GammaAlpha3DSourceTerm::getRethetat(const CFreal Tu, bool prime) 
{
   //if (!(Tu >= 0.0)) CFLog(INFO, "Tu problem: " << Tu << "\n");
  //cf_assert(Tu >= 0.0);   
    
//  const CFreal overTu = 1.0/Tu;
//           
//  if (Tu<=1.3) 
//  {
//    m_Rethetat = (1173.51-589.428*Tu + 0.2196*overTu*overTu);
//  }
//  else 
//  {
//    const CFreal lamco5   = Tu - 0.5658;
//    const CFreal pwtu   = -0.671;
//    
//    m_Rethetat = 331.5*std::pow(lamco5,pwtu);
//  }
  CFreal Rethetat;
  
  if (!prime)
  {
    Rethetat = 0.664*sqrt(400094.0*pow(Tu,-1.38) - 105254.0*pow(Tu,-7.0/8.0));
  }
  else
  {
    Rethetat = 0.332/sqrt(400094.0*pow(Tu,-1.38) - 105254.0*pow(Tu,-7.0/8.0))*(-1.38*400094.0*pow(Tu,-2.38)+105254.0*7.0/8.0*pow(Tu,-15.0/8.0));
  }
  
  return Rethetat;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

CFreal GammaAlpha3DSourceTerm::getRethetatwithPressureGradient(const CFreal avAlpha, const CFreal rhoInfLocal, const CFreal uInfLocal, const CFreal muInfLocal, const CFreal mu, const CFreal rho, const CFreal dpx, const CFreal dpy, const CFreal dpz)
{  
  const CFreal alphaTerm = rhoInfLocal*uInfLocal*uInfLocal/(avAlpha*sqrt(rhoInfLocal*muInfLocal));
  CFreal ReThetat0 = alphaTerm*0.2247141166;
  CFreal ReThetat1 = 0.0;
  const CFreal avV = m_solPhysData[EulerTerm::V]; 
  const CFreal u = m_solPhysData[EulerTerm::VX];
  const CFreal v = m_solPhysData[EulerTerm::VY]; 
  const CFreal w = m_solPhysData[EulerTerm::VZ];
  CFreal dpds = 1.0/avV*(u*dpx + v*dpy + w*dpz);
  const CFreal lambdaTerm1 = -mu/(rho*rho*uInfLocal*uInfLocal*uInfLocal);

  const CFreal dpdsLimit = -m_lambdaLim/(lambdaTerm1*ReThetat0*ReThetat0);
  dpds = min(max(dpds,-dpdsLimit),dpdsLimit);

  const CFreal lambdaTerm = lambdaTerm1*dpds;
  
  const CFuint MAXITER   = 10;
  const CFreal TOL       = 1.0e-6;
  
  for (CFuint iter = 0; iter < MAXITER; ++iter)
  {
    CFreal resReThetat = std::fabs(ReThetat0*TOL);
    
    const CFreal lambda = max(ReThetat0*ReThetat0*lambdaTerm,-0.089);
    const CFreal SLambda = pow(0.09+lambda,0.62);
    //ReThetat1 = alphaTerm*SLambda;

    const CFreal derivDenom = 1.0/pow(0.09+lambda,0.38);

    const CFreal mainF = alphaTerm*SLambda - ReThetat0;
    const CFreal mainFprime = 1.24*alphaTerm*lambdaTerm*ReThetat0*derivDenom - 1.0; 
//CFLog(INFO, "lambdaTerm: " << lambdaTerm << ", alphaTerm: " << alphaTerm << ", Ret0: " << ReThetat0 << ", dpds: " << dpds << ", avV: " << avV << ", lambda: " << lambda << "\n");
    ReThetat1 = max(ReThetat0 - mainF/mainFprime,10.0);
    
    if ( fabs(ReThetat0-ReThetat1) < resReThetat || iter == MAXITER-1) 
    { 
      if (iter == MAXITER-1) CFLog(INFO, "Max iter thetat reached!, Re0: " << ReThetat0 << ", Re1: " << ReThetat1 << ", mainF: " << mainF << ", derivF: " << mainFprime << ", alpha: " << avAlpha << "\n");
          
      break;
    }
    
    if (ReThetat0 != ReThetat0 || ReThetat1 != ReThetat1) 
    { 
      CFLog(INFO, "Nan detected!, Re0: " << ReThetat0 << ", Re1: " << ReThetat1 << ", mainF: " << mainF << ", derivF: " << mainFprime << ", alpha: " << avAlpha << "\n");
          
      break;
    }
    
    ReThetat0 = ReThetat1;
  }
  
  return ReThetat1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

CFreal GammaAlpha3DSourceTerm::getTuInfLocal(const CFreal Rethetat, const CFreal MInfLocal, const CFreal TuLocal)
{  
  const CFreal fM = sqrt(1.0+0.38*pow(MInfLocal,0.6));
  CFreal Tu0 = TuLocal;
  CFreal Tu1 = 0.0;
  
  const CFuint MAXITER   = 50;
  const CFreal TOL       = 1.0e-5;
  
  for (CFuint iter = 0; iter < MAXITER; ++iter)
  {
    CFreal resTu = fabs(Tu0*TOL);
    
    const CFreal mainF = fM*getRethetat(Tu0,false) - Rethetat;
    const CFreal mainFprime = fM*getRethetat(Tu0,true);   
    
    Tu1 = min(max(Tu0 - mainF/mainFprime,0.027),13.0);
    
    //if (!(Tu1>=0.0)) CFLog(INFO, "Tu1: " << Tu1 << ", Tu0: " << Tu0 << ", fM: " << fM << ", mainFPrime: " << mainFprime << "\n");
    
    if ( fabs(Tu0-Tu1) < resTu || iter == MAXITER-1) 
    { 
      if (iter == MAXITER-1) CFLog(INFO, "Max iter TuInfLocal reached! Tu0: " << Tu0 << ", Tu1: " << Tu1 << ", Re: " << Rethetat << ", MInf: " << MInfLocal << ", TuLocal: " << TuLocal << "\n");
          
      break;
    }
    
    Tu0 = Tu1;
  }
  
  return Tu1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void  GammaAlpha3DSourceTerm::getSToStateJacobian(const CFuint iState)
{
  // reset the jacobian
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_stateJacobian[iEq] = 0.0;
  }
  
  SafePtr< NavierStokes3DKLogOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes3DKLogOmega >();
  
  // Set the wall distance before computing the turbulent viscosity
  navierStokesVarSet->setWallDistance(m_currWallDist[iState]);
    
  /// destruction term of k
    
  const CFreal avK     = max((*((*m_cellStates)[iState]))[4],0.0);
  const CFreal avOmega = exp((*((*m_cellStates)[iState]))[5]);
  const CFreal T = max(1.0e-3,(*((*m_cellStates)[iState]))[3]);
  const CFreal p = max(1.0e-3,(*((*m_cellStates)[iState]))[0]);
  const CFreal rho = max(navierStokesVarSet->getDensity(*((*m_cellStates)[iState])),1.0e-8);
  const CFreal R = m_eulerVarSet->getModel()->getR();
  const CFreal overRT = 1.0/(R*T);
  const CFreal pOverRTT = p/(R*T*T);
  const CFreal overOmega = 1.0/avOmega;
  const CFreal avV = sqrt((*((*m_cellStates)[iState]))[1]*(*((*m_cellStates)[iState]))[1]+(*((*m_cellStates)[iState]))[2]*(*((*m_cellStates)[iState]))[2]);
  
  const CFreal avGa    = min(max((*((*m_cellStates)[iState]))[6],0.01),0.99);
    
  const CFreal avAlpha    = max((*((*m_cellStates)[iState]))[7],0.0);
  
  const CFreal mu = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iState]));
  const CFreal mut = navierStokesVarSet->getTurbDynViscosityFromGradientVars(*((*m_cellStates)[iState]), m_cellGrads[iState]);
  
  if (m_isAxisymmetric)
  {
    m_overRadius = 1.0/max(((*m_cellStates)[iState]->getCoordinates())[YY],1.0e-4);
    m_vOverRadius = m_overRadius*(*((*m_cellStates)[iState]))[2];
  }
  
  //Compute Strain 
  getStrain(iState);
  
  // Get Vorticity
  getVorticity(iState);
  
  const CFreal u = (*((*m_cellStates)[iState]))[1];
  const CFreal v = (*((*m_cellStates)[iState]))[2];
  const CFreal dpx = (*(m_cellGrads[iState][0]))[XX];
  const CFreal dpy = (*(m_cellGrads[iState][0]))[YY];
  const CFreal dpz = (*(m_cellGrads[iState][0]))[ZZ];
    
  // compute local freestream values with isentropic assumption
  const CFreal MInf = m_eulerVarSet->getModel()->getMachInf();
  const CFreal pInf = m_eulerVarSet->getModel()->getStaticPressInf();//getPressInf();
  const CFreal gammaIsentropic = m_eulerVarSet->getModel()->getGamma();
  const CFreal uInf = m_eulerVarSet->getModel()->getVelInf();
  const CFreal rhoInf = (MInf*MInf)/(uInf*uInf)*gammaIsentropic*pInf;
  
  const CFreal pTotTerm = 1.0+(gammaIsentropic-1.0)/2.0*MInf*MInf;
  const CFreal pTotExponent = gammaIsentropic/(gammaIsentropic-1.0);
  const CFreal pTotalInf = pow(pTotTerm,pTotExponent)*p;
      
  const CFreal rhoInfLocal = rhoInf*pow(p/pInf,1.0/gammaIsentropic);
  const CFreal MInfLocal = sqrt(2.0/(gammaIsentropic-1.0)*(pow(pTotalInf/p,1.0/pTotExponent)-1.0));
  const CFreal uInfLocal = sqrt(gammaIsentropic*p/rhoInfLocal)*MInfLocal;
  const CFreal TInfLocal = p/(rhoInfLocal*R);
    
  const CFreal Tu = min(max(100.0 * (std::sqrt(2.0*avK/3.0))/(avV),0.027),13.0);
      
  (*((*m_cellStates)[iState]))[3] = TInfLocal;
    
  const CFreal muInfLocal = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iState]));
    
  (*((*m_cellStates)[iState]))[3] = T;
      
  CFreal ReThetat;

  if (!m_PGrad)
  {
    ReThetat = rhoInfLocal*uInfLocal*uInfLocal*0.2247141166/(avAlpha*sqrt(rhoInfLocal*muInfLocal));
  }
  else
  {
    ReThetat = getRethetatwithPressureGradient(avAlpha,rhoInfLocal,uInfLocal,muInfLocal,mu,rho,dpx,dpy,dpz);
  }

  const CFreal TuInfLocal = getTuInfLocal(ReThetat,MInfLocal,Tu);
  
  const CFreal dux = (*(m_cellGrads[iState][1]))[XX];
  const CFreal duy = (*(m_cellGrads[iState][1]))[YY]; 
  const CFreal dvx = (*(m_cellGrads[iState][2]))[XX]; 
  const CFreal dvy = (*(m_cellGrads[iState][2]))[YY]; 
  const CFreal dgammax = (*(m_cellGrads[iState][6]))[XX];
  const CFreal dgammay = (*(m_cellGrads[iState][6]))[YY];
  const CFreal dkx = (*(m_cellGrads[iState][4]))[XX];
  const CFreal dky = (*(m_cellGrads[iState][4]))[YY];
  const CFreal dax = (*(m_cellGrads[iState][7]))[XX];
  const CFreal day = (*(m_cellGrads[iState][7]))[YY];
  
  const CFreal betaStar = navierStokesVarSet->getBetaStar(*((*m_cellStates)[iState]));
  
  const CFreal DkTerm = avOmega * betaStar;
  
  const CFreal Dk = -rho * avK * DkTerm;
  
  const CFreal dmudg = 0.5*(1-2.0*avGa)*(tanh((avGa-0.25)/0.1)+1) + 0.5/0.1*(1-avGa)*avGa/(pow(cosh((avGa-0.25)/0.1),2)) + 1.0;

  //const CFreal limGamma = max(0.01,avGa);
  const CFreal mutGaMod = (avGa+avGa*(1.0-avGa)*0.5*(1.0+tanh((avGa-0.25)/(0.1))));
    
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

  const CFreal coeffTauMu = navierStokesVarSet->getModel().getCoeffTau();
  //const CFreal mutTerm = coeffTauMu*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)-(dux+dvy-m_vOverRadius)*m_vOverRadius)+(duy+dvx)*(duy+dvx));
  const CFreal gammaTerm = avGa + avGa*(1.0-avGa);
  
  CFreal mutTerm;
  CFreal Pk;
  
  if (!m_isSSTV)
  {
    mutTerm = coeffTauMu*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)-(dux+dvy-m_vOverRadius)*m_vOverRadius)+(duy+dvx)*(duy+dvx));

    Pk = (mutTerm*mut - (2./3.)*(avK * rho * gammaTerm)*(dux+dvy+m_vOverRadius));
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
  
  //const CFreal Pk = (mutTerm*mut - (2./3.)*(avK * rho * gammaTerm)*(dux+dvy+m_vOverRadius));
  
  const CFreal twoThirdduxduy = (2./3.)*(dux+dvy+m_vOverRadius);
  
  if (Pk > 0.0)
  {
    CFreal tempSTTerm = 0.0;
       
    if (!m_blockDecoupled)
    {
      //p
      tempSTTerm = mutTerm*avK*overRT*overOmega*mutGaMod;
      if (!m_neglectSSTVTerm) tempSTTerm += -twoThirdduxduy*avK*overRT * gammaTerm;
      
      m_stateJacobian[0][4] += tempSTTerm;
      m_stateJacobian[0][3] -= tempSTTerm;
  
      //T
      tempSTTerm = -mutTerm*avK*overOmega*mutGaMod*pOverRTT;
      if (!m_neglectSSTVTerm) tempSTTerm += twoThirdduxduy*avK*pOverRTT * gammaTerm;
      
      m_stateJacobian[3][4] += tempSTTerm;
      m_stateJacobian[3][3] -= tempSTTerm;
 
      //gamma
      tempSTTerm = mutTerm*mut*dmudg/mutGaMod;
      if (!m_neglectSSTVTerm) tempSTTerm += -twoThirdduxduy * avK * rho*2.0*(1.0-avGa);
      
      m_stateJacobian[6][4] += tempSTTerm;
      m_stateJacobian[6][3] -= tempSTTerm;
      
      if (m_isAxisymmetric && !m_isSSTV)
      {
        //v
        tempSTTerm = mut*coeffTauMu*4./3.*m_overRadius*(-(dux+dvy)+v*2.0)- (2./3.)*(avK * rho * gammaTerm)*m_overRadius;
        m_stateJacobian[2][4] += tempSTTerm;
        m_stateJacobian[2][3] -= tempSTTerm;
      }
    }
  
    //k
    tempSTTerm = mutTerm*rho*overOmega*mutGaMod;
    if (!m_neglectSSTVTerm) tempSTTerm += -twoThirdduxduy*rho * gammaTerm;
    
    m_stateJacobian[4][4] += tempSTTerm;
    m_stateJacobian[4][3] -= tempSTTerm;
  
    //logOmega
    tempSTTerm = -mutTerm*rho*avK*overOmega*mutGaMod;
    m_stateJacobian[5][4] += tempSTTerm;
    m_stateJacobian[5][3] -= tempSTTerm;
  }
  
  /// production of logOmega
  
  const CFreal gamma = navierStokesVarSet->getGammaCoef();
  const CFreal blendingCoefF1 = navierStokesVarSet->getBlendingCoefficientF1();
  const CFreal sigmaOmega2 = navierStokesVarSet->getSigmaOmega2();
  
  const CFreal pOmegaFactor = (1. - blendingCoefF1) * 2. * sigmaOmega2 * ((*(m_cellGrads[iState][4]))[XX]*(*(m_cellGrads[iState][5]))[XX] + (*(m_cellGrads[iState][4]))[YY]*(*(m_cellGrads[iState][5]))[YY]);
  
  const CFreal Pomega = (gamma*rho/mut) * Pk * overOmega + rho * overOmega * pOmegaFactor;
  
  if (Pomega > 0.0)
  {
    if (!m_blockDecoupled)
    {
      //p
      m_stateJacobian[0][5] += gamma*(mutTerm*overRT*overOmega - twoThirdduxduy*overRT/mutGaMod*gammaTerm) + pOmegaFactor*overRT*overOmega;
  
      //T
      m_stateJacobian[3][5] += gamma*pOverRTT*(-mutTerm*overOmega + twoThirdduxduy/mutGaMod*gammaTerm) - pOmegaFactor*pOverRTT*overOmega;
    }
 
    //logOmega
    m_stateJacobian[5][5] +=  -gamma*rho*mutTerm*overOmega - pOmegaFactor*rho*overOmega;
  }
  
  // P gamma
  const CFreal cfg1 = 1.735;
  const CFreal cfg2 = 5.45;
  const CFreal cfg3 = 0.95375;
  const CFreal cfg4 = 2.2;
    
  CFreal fg = 1.0;
    
  if (avGa < 0.45) fg -= exp(-cfg1*tan(cfg2*avGa-cfg3)-cfg4);
    
  const CFreal fMnsigmaTerm = 1.0+0.58*pow(MInfLocal,0.6);
  const CFreal fMnsigma = 1.0/(fMnsigmaTerm*fMnsigmaTerm);
    
  CFreal dpds = 1.0/avV*(u*dpx + v*dpy);

  CFreal PRC = 1.0;
    
  if (m_PGrad)
  {
    const CFreal lambdaTerm1 = -mu/(rho*rho*uInfLocal*uInfLocal*uInfLocal);

    const CFreal dpdsLimit = -m_lambdaLim/(lambdaTerm1*ReThetat*ReThetat);
    dpds = min(max(dpds,-dpdsLimit),dpdsLimit);
    
    const CFreal kPGrad = -muInfLocal/(rhoInfLocal*rhoInfLocal*uInfLocal*uInfLocal*uInfLocal)*fabs(1.0-MInfLocal*MInfLocal)*dpds;
    
      
    if (kPGrad < 0)
    {
      const CFreal PRCfactor1 = 474.0*pow(TuInfLocal,-2.9);
      const CFreal PRCexponent = 1.0-exp(2.0e6*kPGrad);
      PRC = pow(PRCfactor1,PRCexponent);
    }
    else
    {
      const CFreal PRCexponent = -3227.0*pow(kPGrad,0.5985);
      PRC = pow(10.0,PRCexponent);
    }
  }
    
  const CFreal fk = 1.0+fg*(PRC-1.0); 
   
  const CFreal nsigma = 1.25e-11*pow(TuInfLocal,7.0/4.0)*fk*fMnsigma;
  
  const CFreal betaGA = sqrt(nsigma)*uInfLocal*rhoInfLocal/muInfLocal;
  
  const CFreal gTerm = (1-avGa)*sqrt(-log(1.0-avGa));
  
  const CFreal brc2 = betaGA*rho*avV*2.0;
  
  CFreal prodTerm_Ga = fg*gTerm*brc2;
  
  if (prodTerm_Ga > 0.0)
  {    
    const CFreal dgTermdg = (2.0*log(1-avGa)+1.0)/(2.0*sqrt(-log(1-avGa))); 
    
    CFreal dfgdg = 0.0;
    
    if (avGa<0.45) dfgdg = cfg1*cfg2/pow(cos(cfg2*avGa-cfg3),2) * exp(cfg1*tan(cfg3-cfg2*avGa)-cfg4);
    
    //const CFreal dbdg = 0.5*uInfLocal*rhoInfLocal/muInfLocal/sqrt(nsigma) * 1.25e-11*pow(TuInfLocal,7.0/4.0) * fMnsigma * (PRC-1.0) * dfgdg;
     
    //gamma
    //m_stateJacobian[6][6] += brc2 * fg * dgTermdg + brc2 * gTerm * dfgdg;// + gTerm * fg * 2.0 * rho * avV * dbdg;
   
    if (!m_blockDecoupled)
    {
      const CFreal gTermFgBeta2 = gTerm * fg * betaGA *2.0;
      
      //p
      m_stateJacobian[0][6] +=  gTermFgBeta2 * overRT * avV;
  
      //T
      m_stateJacobian[3][6] +=  -gTermFgBeta2 * pOverRTT * avV;
      
      //u
      m_stateJacobian[1][6] +=  gTermFgBeta2 * rho * u/avV;
              
      //v
      m_stateJacobian[2][6] +=  gTermFgBeta2 * rho * v/avV;    
    }
  }
  
  // D gamma
  const CFreal cEg = 20.0/0.57;
  
  const CFreal fMuGamma = 1.0-exp(-256.0*(m_currWallDist[iState]*uInfLocal*rhoInfLocal/muInfLocal)*(m_currWallDist[iState]*uInfLocal*rhoInfLocal/muInfLocal));
  const CFreal fMMuGamma = (1.0+0.26*(gammaIsentropic-1.0)/2.0*MInfLocal*MInfLocal)*sqrt(1+0.38*pow(MInfLocal,0.6));
  
  //const CFreal gammaLim = std::min(std::max(0.01,avGa),0.99);
  const CFreal muGamma = 0.57*pow(-log(1.0-avGa),-5.0/6.0*(1.0-avGa))*fMuGamma*fMMuGamma*mu;
    
  const CFreal dudn = -1.0/(avV*avV)*(u*u*duy-v*v*dvx+u*v*(dvy-dux));//-duy; //
  const CFreal dgammadn = 1.0/avV*(v*dgammax - u*dgammay);//-dgammay; //
  
  const CFreal dgTerm = -cEg * avV/(uInfLocal*uInfLocal) * dudn * dgammadn;
    
  CFreal  destructionTerm_Ga  = dgTerm * muGamma;
    
  //destructionTerm_Ga = max(-10.0*fabs(prodTerm_Ga),destructionTerm_Ga);
  //if (avGa<0.4) destructionTerm_Ga = max(-fabs(prodTerm_Ga),destructionTerm_Ga);
  
  if (destructionTerm_Ga < 0.0 && m_addDGDA)
  {
    const CFreal mgdg = muGamma * 5.0/6.0 * (log(-log(1-avGa)) + 1.0/log(1-avGa));
    
    // gamma
    //m_stateJacobian[6][6] +=  dgTerm * mgdg;
  
    if (!m_blockDecoupled)
    {
      // u
      m_stateJacobian[1][6] +=  destructionTerm_Ga * u/(avV*avV);
      
      // v
      m_stateJacobian[2][6] +=  destructionTerm_Ga * v/(avV*avV);      
    }
  }
  
  // production  Alpha
  
  const CFreal cpa1 = 0.03;
  //const CFreal cpa2 = 50.0;

  CFreal Rethetac = getRethetat(Tu, false);
    
  ///@todo check if this is local M or local MInf
  Rethetac *= sqrt(1.0+0.38*pow(MInfLocal,0.6));
    
  CFreal lambda = 0.0;

  if (m_PGrad)
  {
    lambda = -Rethetac*Rethetac*mu/(rho*rho*uInfLocal*uInfLocal*uInfLocal)*dpds;
  }
    
  const CFreal alphac = rhoInfLocal*uInfLocal*uInfLocal*pow(0.09+lambda,0.62)/(Rethetac*sqrt(rhoInfLocal*muInfLocal));
    
  const CFreal t = (500.0 * mu )/(rho * avV * avV);
  
  const CFreal dkdn = 1.0/avV*(v*dkx - u*dky);//-dky; //  
    
  const CFreal  Rew         = (rho * m_currWallDist[iState] * m_currWallDist[iState] * avOmega)/(mu);   
  const CFreal  Fwake1      = (1.0e-5 * Rew)*(1.0e-5 * Rew);   
  const CFreal  Fwake       = exp(-Fwake1);
  const CFreal  thetaBL     = (ReThetat*mu)/(rho*avV);
     
  const CFreal  delta       = (375.0 * m_vorticity * m_currWallDist[iState] * thetaBL)/(avV);
    
  const CFreal  coefFtheta0 = (m_currWallDist[iState]/delta)*(m_currWallDist[iState]/delta)*(m_currWallDist[iState]/delta)*(m_currWallDist[iState]/delta);
  const CFreal  coefFtheta1 = exp(-coefFtheta0);
  const CFreal  Ftheta1     = Fwake * coefFtheta1;
  const CFreal  Ftheta3     = 1.0-(((avGa-0.01)/(0.99-0.01))*((avGa-0.01)/(0.99-0.01)));
  const CFreal  Ftheta4     = std::max(Ftheta1,Ftheta3);
  CFreal  Fthetat     = std::min(Ftheta4,1.0);
  
  const CFreal dkdnFactor = fabs(dkdn)/max(avK,1.0e-6);

  if (dkdnFactor>100.0) Fthetat = 1.0;
    
  CFreal prodTerm_alpha = cpa1 * (rho/t) * (alphac - avAlpha) * (1.0 - Fthetat);

  // compute destruction term of alpha
  const CFreal cna = 0.4;
  const CFreal dadn = 1.0/avV*(v*dax - u*day);
    
  CFreal destructionTerm_alpha = -cna*rho*avV*dadn;
  
  ///@todo Fthetat and alpha_c comtribution
  if (prodTerm_alpha > 0.0 || !m_limPRe)
  {
    // alpha
    m_stateJacobian[7][7] +=  -cpa1 * (rho/t) * (1.0 - Fthetat);
  
    if (!m_blockDecoupled)
    {
      const CFreal tTerm = 2.0*rho*avV*avV/(500.0*mu);
        
      // p
      m_stateJacobian[0][7] +=  cpa1 * tTerm * overRT * (alphac - avAlpha) * (1.0 - Fthetat);
  
      // T
      m_stateJacobian[3][7] +=  -cpa1 * tTerm * pOverRTT * (alphac - avAlpha) * (1.0 - Fthetat);
  
      // u
      m_stateJacobian[1][7] +=  cpa1 * 2.0*rho*rho*u/(500.0*mu) * (alphac - avAlpha) * (1.0 - Fthetat);
  
      // v
      m_stateJacobian[2][7] +=  cpa1 * 2.0*rho*rho*v/(500.0*mu) * (alphac - avAlpha) * (1.0 - Fthetat);
    }
  }
  
  if (destructionTerm_alpha < 0.0 && !m_blockDecoupled)
  {
      const CFreal daTerm = -cna*dadn;
        
      // p
      m_stateJacobian[0][7] +=  daTerm*avV*overRT;
  
      // T
      m_stateJacobian[3][7] +=  -daTerm*avV*pOverRTT;
  
      // u
      m_stateJacobian[1][7] +=  daTerm*rho*u/avV;
  
      // v
      m_stateJacobian[2][7] +=  daTerm*rho*v/avV;
  }
  
  if (m_isAxisymmetric)
  { 
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
      
      //rho gamma
      if (!m_blockDecoupled)
      {
        m_stateJacobian[0][6] += -vrOverRT*avGa;
        m_stateJacobian[2][6] += -rhor*avGa;
        m_stateJacobian[3][6] += vrPOverRTT*avGa;
      }
      m_stateJacobian[6][6] += -rhovr;
      
      //rho alpha
      if (!m_blockDecoupled)
      {
        m_stateJacobian[0][7] += -vrOverRT*avAlpha;
        m_stateJacobian[2][7] += -rhor*avAlpha;
        m_stateJacobian[3][7] += vrPOverRTT*avAlpha;
      }
      m_stateJacobian[7][7] += -rhovr;
  }

  if (m_addUpdateCoeff)
  {
    DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
    CFreal update = 0.0;
  
    const CFreal PkUpdate = max(Pk,0.0);
  
    const CFreal DkUpdate = min(Dk,-1.0e-10);
  
    const CFreal G = sqrt(-PkUpdate*rho/DkUpdate);
  
    const CFreal Gk = pow(G,1.0/gamma);
  
    const CFreal kUpdateTerm = (PkUpdate + DkUpdate)/max(avK,0.001)/(Gk-1) + m_stateJacobian[4][4];
  
    update = max(0.0,kUpdateTerm);
  
    //if (update>1000) CFLog(INFO,"kUpdate: " << update << "\n");
  
    const CFreal logOmega = (*((*m_cellStates)[iState]))[5];
  
    const CFreal Gomega = pow(G,1.0/logOmega);
  
    const CFreal PomegaUpdate = max(Pomega,0.0);
  
    const CFreal DomegaUpdate = min(Domega,0.0);
  
    const CFreal omegaUpdateTerm = (PomegaUpdate + DomegaUpdate)/max(logOmega,0.1)/(Gomega-1) + m_stateJacobian[5][5];
  
    update = max(update,omegaUpdateTerm);
  
    // get the local ID of the current sol pnt
    const CFuint solID = (*m_cellStates)[iState]->getLocalID();
  
    //if (update>1000) CFLog(INFO,"old: " << updateCoeff[solID] << ", new: " << update << ", J: " << m_solPntJacobDets[iState] << "\n\n");
  
    update *= (2.0*m_order+1)*m_solPntJacobDets[iState];//*m_solPntJacobDets[iState];
  
    if (update>0.5*updateCoeff[solID]) CFLog(INFO,"Large influence ST Dt: oldCoeff: " << updateCoeff[solID] << ", ST addition: " << update << ", J: " << m_solPntJacobDets[iState] << "\n\n");
    
    // add the wave speed update previously computed
    updateCoeff[solID] += update;
  }
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
