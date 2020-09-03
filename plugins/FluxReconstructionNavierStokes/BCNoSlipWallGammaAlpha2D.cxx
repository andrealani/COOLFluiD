#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/BCNoSlipWallGammaAlpha2D.hh"

#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCNoSlipWallGammaAlpha2D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNavierStokesModule >
  BCNoSlipWallGammaAlpha2DProvider("BCNoSlipWallGammaAlpha2D");

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallGammaAlpha2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("yWallVelocity","Y-component of a velocity vector of the wall.");
  options.addConfigOption< CFreal >("xWallVelocity","X-component of a velocity vector of the wall.");
  options.addConfigOption< CFreal >("KWall","Wall value for turbulent intensity");
  options.addConfigOption< CFreal,Config::DynamicOption<> >("T","wall static temperature");
  options.addConfigOption< CFreal >("q","wall heat flux");
  options.addConfigOption< bool >("HeatFlux","bool to tell if the wall has constant heat flux (possibly zero), default true.");
  options.addConfigOption< CFuint,Config::DynamicOption<> >("ChangeToIsoT","Iteration after which to switch to an isothermal BC.");
  options.addConfigOption< CFreal >("WallDist","Characteristic distance of first sol pnt from the wall.");
  options.addConfigOption< CFreal >("OmegaWallFactor","Factor by which to multiply omegaWall each iteration until it is the theoretical value (Default 1.01).");
  options.addConfigOption< CFuint >("ImposeOmegaWallIter","Iteration at which to impose theoretical omegaWall value.");
  options.addConfigOption< CFreal,Config::DynamicOption<>  >("GammaThreshold","Threshold beyond which gamma is considered to be transitioning if transition already occurred (Default 0.02).");
  options.addConfigOption< bool,Config::DynamicOption<>  >("EnforceLaminar","Flag to enforce laminar flow at all times (Default false).");
}

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallGammaAlpha2D::BCNoSlipWallGammaAlpha2D(const std::string& name) :
  BCStateComputer(name),
  m_varSetTurb(CFNULL),
  m_diffVarTurb(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_prevIter(0),
  m_transitionOccurred(false)
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);
   
  m_wallT = 0.0;
   setParameter("T",&m_wallT);

  m_wallQ = 0.0;
   setParameter("q",&m_wallQ);

  m_heatFlux= true;
   setParameter("HeatFlux",&m_heatFlux);
   
  m_changeToIsoT = MathTools::MathConsts::CFuintMax();
   setParameter("ChangeToIsoT",&m_changeToIsoT);
   
  m_xWallVelocity = 0.0;
   setParameter("xWallVelocity",&m_xWallVelocity);
  
  m_yWallVelocity = 0.0;
   setParameter("yWallVelocity",&m_yWallVelocity);
  
  m_wallK = 0.0;
   setParameter("KWall",&m_wallK);
   
   m_wallDist = 1.0e-5;
   setParameter("WallDist",&m_wallDist);
   
   m_imposeOmegaWallIter = 0;
   setParameter("ImposeOmegaWallIter",&m_imposeOmegaWallIter);
   
   m_omegaWallFactor = 1.01;
   setParameter("OmegaWallFactor",&m_omegaWallFactor);
   
   m_gammaThreshold = 0.02;
   setParameter("GammaThreshold",&m_gammaThreshold);
   
   m_enforceLaminar = false;
   setParameter("EnforceLaminar",&m_enforceLaminar);
}

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallGammaAlpha2D::~BCNoSlipWallGammaAlpha2D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallGammaAlpha2D::computeGhostStates(const vector< State* >& intStates,
                                                    vector< State* >& ghostStates,
                                                    const std::vector< RealVector >& normals,
                                                    const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());
  
  const CFuint iK = m_varSetTurb->getModel()->getFirstScalarVar(0);
  const CFuint nbTurbVars = m_varSetTurb->getModel()->getNbScalarVars(0);
  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  if (iter >= m_changeToIsoT && m_heatFlux)
  {
    m_heatFlux = false;
  }

  // get some data from the physical model
  const CFreal gamma = m_varSetTurb->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
  const CFreal R = m_varSetTurb->getModel()->getR();

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);

    // set the physical data starting from the inner state
    m_varSetTurb->computePhysicalData(intState,m_intSolPhysData);
    
    if (m_heatFlux)
    {
      // set the physical data for the ghost state
      m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::RHO];
      m_ghostSolPhysData[EulerTerm::VX]  = m_xWallVelocity;//2.0*m_xWallVelocity-m_intSolPhysData[EulerTerm::VX];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::VY]  = m_yWallVelocity;//2.0*m_yWallVelocity-m_intSolPhysData[EulerTerm::VY];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::V] = sqrt(m_ghostSolPhysData[EulerTerm::VX]*m_ghostSolPhysData[EulerTerm::VX]+m_ghostSolPhysData[EulerTerm::VY]*m_ghostSolPhysData[EulerTerm::VY]);
      m_ghostSolPhysData[EulerTerm::P]   = m_intSolPhysData[EulerTerm::P];
      m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*m_ghostSolPhysData[EulerTerm::P]
                                            + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
                                                  m_ghostSolPhysData[EulerTerm::V]*
                                                  m_ghostSolPhysData[EulerTerm::V]
                                         )/m_ghostSolPhysData[EulerTerm::RHO];
      m_ghostSolPhysData[EulerTerm::A] = sqrt(gamma*m_ghostSolPhysData[EulerTerm::P]/m_ghostSolPhysData[EulerTerm::RHO]);
      m_ghostSolPhysData[EulerTerm::T] = m_intSolPhysData[EulerTerm::T];
      m_ghostSolPhysData[EulerTerm::E] = m_ghostSolPhysData[EulerTerm::H] -
                                         (m_ghostSolPhysData[EulerTerm::P]/m_ghostSolPhysData[EulerTerm::RHO]);
    }
    else
    {
      const CFreal innerT = m_intSolPhysData[EulerTerm::P]/(R*m_intSolPhysData[EulerTerm::RHO]);
      
      CFreal ghostT = m_wallT;

      // set the physical data for the ghost state
      m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::P]/(R*ghostT); //m_intSolPhysData[EulerTerm::RHO];
      m_ghostSolPhysData[EulerTerm::VX]  = m_xWallVelocity;//2.0*m_xWallVelocity-m_intSolPhysData[EulerTerm::VX];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::VY]  = m_yWallVelocity;//2.0*m_yWallVelocity-m_intSolPhysData[EulerTerm::VY];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::V] = sqrt(m_ghostSolPhysData[EulerTerm::VX]*m_ghostSolPhysData[EulerTerm::VX]+m_ghostSolPhysData[EulerTerm::VY]*m_ghostSolPhysData[EulerTerm::VY]);
      m_ghostSolPhysData[EulerTerm::P]   = m_intSolPhysData[EulerTerm::P]; //ghostP;
      m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*m_ghostSolPhysData[EulerTerm::P]
                                            + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
                                                  m_ghostSolPhysData[EulerTerm::V]*
                                                  m_ghostSolPhysData[EulerTerm::V]
                                         )/m_ghostSolPhysData[EulerTerm::RHO];
      m_ghostSolPhysData[EulerTerm::A] = sqrt(gamma*m_ghostSolPhysData[EulerTerm::P]/m_ghostSolPhysData[EulerTerm::RHO]);
      m_ghostSolPhysData[EulerTerm::T] = ghostT;
      m_ghostSolPhysData[EulerTerm::E] = m_ghostSolPhysData[EulerTerm::H] -
                                         (m_ghostSolPhysData[EulerTerm::P]/m_ghostSolPhysData[EulerTerm::RHO]);
    }

    m_ghostSolPhysData[iK] = m_wallK;
    
    const CFreal pdim =  m_ghostSolPhysData[EulerTerm::P] * m_varSetTurb->getModel()->getPressRef();
    const CFreal Tdim =  m_ghostSolPhysData[EulerTerm::T] * m_varSetTurb->getModel()->getTempRef();
    const CFreal mu = m_diffVarTurb->getModel().getDynViscosityDim(pdim, Tdim)/(m_diffVarTurb->getModel().getReferencePhysicalData())[NSTurbTerm::MU];
    
    // check if it is k-omega and not SA
    if(nbTurbVars == 2 || nbTurbVars == 4)
    {
      //Compute distance to innerstate
      CFreal y0 = m_wallDist;//1.e-9;
    
      //avoid too small distances
      //y0 = std::max(y0, 10.e-10);
    
      CFreal nu = mu / m_ghostSolPhysData[EulerTerm::RHO];
    
      //this is not the best, but it avoids having to code another BC! because I
      //would have to dynamic cast to the KOmega varset to get the beta1
      const CFreal beta1 = 0.075;
      
      const CFuint omegaID = iK+1; 
      
      ///@todo here should this be adimensionalized (by the distance)???
      //Menter's definition
      // for stability, gradually increase w_wall
      const CFreal omegaWallTh = log((60. * nu) / (beta1 * y0 * y0));//(10. * 6. * nu) / (beta1 * y0 * y0);
      //const CFreal omegaWall = iter < m_imposeOmegaWallIter ? min(omegaWallTh,m_omegaWallFactor*m_intSolPhysData[omegaID]): omegaWallTh;
      const CFreal omegaWall = iter < m_imposeOmegaWallIter ? min(omegaWallTh,log(m_omegaWallFactor)*iter+5.0): omegaWallTh;

      if (m_prevIter < iter && iter < m_imposeOmegaWallIter) 
      {
        CFLog(INFO, "OmegaWall log difference (-inf -> 0): " << log10(fabs(1+(exp(omegaWall)-exp(omegaWallTh))/exp(omegaWallTh))) << "\n");
        
        m_prevIter = iter;
      }
      
      cf_assert(omegaWall>0.0);

      m_ghostSolPhysData[omegaID] = omegaWall;
      
      m_ghostSolPhysData[EulerTerm::E] += m_ghostSolPhysData[iK];
      m_ghostSolPhysData[EulerTerm::H] += m_ghostSolPhysData[iK];
    }
    
    // check if LCTM is active
    if (nbTurbVars == 4)
    {
//       const CFreal muTot = m_diffVarTurb->getDynViscosity(intStates, m_cellGradFlxPnt[flxIdx]);
//        const CFreal tau = m_muWall*(m_unitNormalFlxPnts[flxIdx][YY]*MathFunctions::innerProd(*(m_cellGradFlxPnt[flxIdx][m_UID]),m_unitNormalFlxPnts[flxIdx]) -
//       m_unitNormalFlxPnts[flxIdx][XX]*MathFunctions::innerProd(*(m_cellGradFlxPnt[flxIdx][m_VID]),m_unitNormalFlxPnts[flxIdx]))
//        const CFreal wallShear = 
//        
//        m_muWall = m_diffVar->getDynViscosity(*(m_cellStatesFlxPnt[flxIdx]), m_cellGradFlxPnt[flxIdx]);
//  
//  // this will not work adimensional
//  // const CFreal tauRef = 1.0; // mu_ref*Uref/refLength;
//  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
//  const bool adim = Framework::PhysicalModelStack::getActive()->getImplementor()->isAdimensional();
//  const CFuint dim = PhysicalModelStack::getActive()->getDim();
//  
//  m_tau = 0.0;
//  m_Cf = 0.0;
//  
//  if (dim == DIM_2D) {
//    
//    m_tau = m_muWall*
//      (m_unitNormalFlxPnts[flxIdx][YY]*MathFunctions::innerProd(*(m_cellGradFlxPnt[flxIdx][m_UID]),m_unitNormalFlxPnts[flxIdx]) -
//       m_unitNormalFlxPnts[flxIdx][XX]*MathFunctions::innerProd(*(m_cellGradFlxPnt[flxIdx][m_VID]),m_unitNormalFlxPnts[flxIdx]));
//    
//    if(adim){
//      //check this
//      const CFreal mInf = m_updateVarSet->getModel()->getMachInf();
//      const CFreal Re = m_diffVar->getModel().getReynolds();
//      m_Cf = m_tau / (0.5*mInf*sqrt(gamma)*Re);
//    }
//    else{
//      // friction coefficient 
//      m_Cf = m_tau / (0.5*m_rhoInf*m_uInf*m_uInf);
//    } 
//    
//    m_frictionForces[XX] =  m_Cf*m_unitNormalFlxPnts[flxIdx][YY];
//    m_frictionForces[YY] = -m_Cf*m_unitNormalFlxPnts[flxIdx][XX];
//  }
        
      m_ghostSolPhysData[EulerTerm::GAMMA] = m_intSolPhysData[EulerTerm::GAMMA];
      //CFLog(INFO, "is turb: " << m_transitionCriterion[iState] << ", transition occurred: " << m_transitionOccurred << "\n");
      if ((m_transitionCriterion[iState] || (m_intSolPhysData[iK+2] > m_gammaThreshold)) && !m_enforceLaminar)// && m_transitionOccurred))
      {
        // gamma
        m_ghostSolPhysData[iK+2] = m_intSolPhysData[iK+2];
        
        if (!m_transitionOccurred)
        {
          m_transitionOccurred = true;
          CFLog(INFO,"\n\n\n Transition occurred!\n\n\n");
        }
      }
      else
      {
        // gamma
        m_ghostSolPhysData[iK+2] = 0.01;
      }
      
      // Alpha
      m_ghostSolPhysData[iK+3] = m_intSolPhysData[iK+3];
    }
    else if (nbTurbVars == 3)
    {
      m_ghostSolPhysData[EulerTerm::GAMMA] = m_intSolPhysData[EulerTerm::GAMMA];

      // gamma
      m_ghostSolPhysData[iK+1] = m_intSolPhysData[iK+1];
      
      // Re
      m_ghostSolPhysData[iK+2] = m_intSolPhysData[iK+2]; 
    }

    // set the ghost state from its physical data
    m_varSetTurb->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
    
    //CFLog(INFO, "ghostState: " << ghostState << "\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallGammaAlpha2D::computeGhostGradients
                                                    (const std::vector< std::vector< RealVector* > >& intGrads,
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
  
  if (m_heatFlux)
  { 
    for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
    {
      // normal
      const RealVector& normal = normals[iState];
    
      // temperature
      RealVector& tempGradI = *intGrads  [iState][3];
      RealVector& tempGradG = *ghostGrads[iState][3];
      const CFreal nTempGrad = tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY];
      tempGradG = tempGradI - nTempGrad*normal + m_wallQ*normal; //tempGradI - 2.0*nTempGrad*normal + m_wallQ*normal;
      
      // pressure
//      RealVector& pGradI = *intGrads  [iState][0];
//      RealVector& pGradG = *ghostGrads[iState][0];
//      const CFreal nPGrad = pGradI[XX]*normal[XX] + pGradI[YY]*normal[YY];
//      pGradG = pGradI - nPGrad*normal;
    }
  }
  
  // check if LCTM is active
  if (nbrGradVars > 6)
  {
    for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
    {
      // normal
      const RealVector& normal = normals[iState];
    
      // gamma
      RealVector& gammaGradI = *intGrads  [iState][nbrGradVars-2];
      RealVector& gammaGradG = *ghostGrads[iState][nbrGradVars-2];
      const CFreal nGammaGrad = gammaGradI[XX]*normal[XX] + gammaGradI[YY]*normal[YY];
      gammaGradG = gammaGradI - nGammaGrad*normal; //tempGradI - 2.0*nTempGrad*normal + m_wallQ*normal;
      
      // Ret
      RealVector& RetGradI = *intGrads  [iState][nbrGradVars-1];
      RealVector& RetGradG = *ghostGrads[iState][nbrGradVars-1];
      const CFreal nRetGrad = RetGradI[XX]*normal[XX] + RetGradI[YY]*normal[YY];
      RetGradG = RetGradI - nRetGrad*normal;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallGammaAlpha2D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  m_varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb2DVarSet>();

  m_diffVarTurb = getMethodData().getDiffusiveVar().d_castTo<DiffTurb2DVarSet>();

  m_varSetTurb->getModel()->resizePhysicalData(m_intSolPhysData);
  m_varSetTurb->getModel()->resizePhysicalData(m_ghostSolPhysData);
  
  m_xWallVelocity /= m_varSetTurb->getModel()->getVelRef();
  m_yWallVelocity /= m_varSetTurb->getModel()->getVelRef();
  
  cf_assert(m_wallK >= 0.0);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

