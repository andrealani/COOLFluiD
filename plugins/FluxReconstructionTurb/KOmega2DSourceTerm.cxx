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

CFreal KOmega2DSourceTerm::GetNSSourceTerm()  
{ 
//  using namespace COOLFluiD::Framework;
//  using namespace COOLFluiD::Physics::NavierStokes;
//    
//  const CFuint uID = getStateVelocityIDs()[XX];
//  const CFuint vID = getStateVelocityIDs()[YY];
//  const CFreal coeffTauMu = _diffVarSet->getModel().getCoeffTau();
//  const CFreal Tau_tt = (-2./3.)*coeffTauMu*((*(_gradients[uID]))[XX] + (*(_gradients[vID]))[YY] - 2*_vOverRadius);
//  const CFreal Source3 = _physicalData[EulerTerm::P] - Tau_tt;
//  return Source3;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::prepareComputeSource()
{  
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Physics::NavierStokes;
  
  //DataHandle<CFreal> normals  = this->socket_normals.getDataHandle();
  //DataHandle<CFint> isOutward =  this->socket_isOutward.getDataHandle();
  //DataHandle<CFreal> volumes  = socket_volumes.getDataHandle();
  
  cf_assert(m_eulerVarSet.isNotNull());
  
  //m_isAxisymmetric = this->getMethodData().isAxisymmetric();
  
//  _Radius = (_isAxisymmetric) ? (currState->getCoordinates())[YY] : 1.; 
//  if((_Radius > MathTools::MathConsts::CFrealEps()) && (_isAxisymmetric)) {
//    _vOverRadius = _physicalData[EulerTerm::VX]/_Radius;
//  }
}
      
////////////////////////////////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::addSourceTerm(RealVector& resUpdates)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Physics::NavierStokes; 
  
//  for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
//  {
//    m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]), m_solPhysData);
//    
//    const CFreal R = m_eulerVarSet->getModel()->getR();
//
//    CFreal dUdX = 0.0;
//    CFreal dVdR = 0.0;
//    
//    const CFuint uID = 1;
//    const CFuint vID = 2;
//    
//    if (true) //Puvt
//    {
//      dUdX = (*(m_cellGrads[iSol]))[uID][XX];
//      dVdR = (*(m_cellGrads[iSol]))[vID][YY];
//    }
//    else
//    {
//      const CFreal invRho = 1.0/(*((*m_cellStates)[iSol]))[0];
//      const CFreal u = invRho*(*((*m_cellStates)[iSol]))[uID];
//      const CFreal v = invRho*(*((*m_cellStates)[iSol]))[vID];
//      
//      // apply chain rule
//      dUdX = invRho*((*(m_cellGrads[iSol]))[uID][XX] - u*(*(m_cellGrads[iSol]))[0][XX]);
//      dVdR = invRho*((*(m_cellGrads[iSol]))[vID][YY] - v*(*(m_cellGrads[iSol]))[0][YY]);
//    }
//    
//    const CFreal avV = m_solPhysData[EulerTerm::VY];
//    
//    // @todo this will not work if gradients are needed (Menter SST turb model)
//    const CFreal mu = navierStokesVarSet->getDynViscosity(*((*m_cellStates)[iSol]), m_dummyGradients);
//    const CFreal coeffMu = navierStokesVarSet->getModel().getCoeffTau()*2.0/3.0*mu;
//    const CFreal invR = 1.0/(m_cell->computeCoordFromMappedCoord((*m_solPntsLocalCoords)[iSol]))[YY];
//    const CFreal tauThetaTheta = -coeffMu*(dUdX + dVdR - 2.0*avV*invR);
//    
//    // AL: check this in Hontzatko's report (dp!)
//    //m_srcTerm[vID] = m_solPhysData[EulerTerm::P] - tauThetaTheta;
//    resUpdates[m_nbrEqs*iSol + vID] = m_solPhysData[EulerTerm::P] - tauThetaTheta;
//
////     for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resID)
////     {
////       rhs[resID] += resFactor*m_solPntJacobDets[iSol]*m_srcTerm[iEq];
////     }
//  }
  
  //prepareComputeSource(); 
  
  // compute PUVTKOmega by averaging the nodes
  // NO!!! If we do it like that we nearly certainly
  // get negative values!!!
  // So we just take the state value
//  const State& avState = *element->getState(0);
  
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
    
//  const CFreal mut = _diffVarSet->getTurbDynViscosityFromGradientVars(avState, _gradients);
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
//  
//  //Computation of the source term
//  const CFuint vID = getStateVelocityIDs()[YY];
//  source[vID] = (_isAxisymmetric) ? GetNSSourceTerm() : 0.0 ;
//  
//  //What we do with the source term depends if
//  //we are computing the jacobian or not
//  const bool isPerturb = this->getMethodData().isPerturb();
//  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();
//  if(isPerturb)
//  {
//    /// Compute the jacobian contribution
//    // only perturb the negative part of the source term
//    if(iPerturbVar == kID)
//    {
//      source[kID] = _destructionTerm_k;
//      source[kID] += _unperturbedPositivePart[0];
//    }
//    else
//    {
//      source[kID] = _unperturbedNegativePart[0];
//      source[kID] += _unperturbedPositivePart[0];
//    }
//
//    if(iPerturbVar == omegaID)
//    {
//      source[omegaID] = _destructionTerm_Omega;
//      source[omegaID] += _unperturbedPositivePart[1];
//    }
//    else
//    {
//      source[omegaID] = _unperturbedNegativePart[1];
//      source[omegaID] += _unperturbedPositivePart[1];
//    }
//  }
//  else
//  {
//    /// Compute the rhs contribution
//    // and Store the unperturbed source terms
//    source[kID] = _prodTerm_k;
//    source[kID] += _destructionTerm_k;
//    _unperturbedPositivePart[0] = _prodTerm_k;
//    _unperturbedNegativePart[0] = _destructionTerm_k;
//
//    source[omegaID] = _prodTerm_Omega;
//    source[omegaID] += _destructionTerm_Omega;
//    _unperturbedPositivePart[1] = _prodTerm_Omega;
//    _unperturbedNegativePart[1] = _destructionTerm_Omega;
//  }
//  
//  // Finally multiply by the cell volume
//  source *= _volumes_elemID;
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
