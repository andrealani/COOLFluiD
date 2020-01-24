#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/BCNoSlipWallTurb2D.hh"

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
    BCNoSlipWallTurb2D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNavierStokesModule >
  BCNoSlipWallTurb2DProvider("BCNoSlipWallTurb2D");

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallTurb2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("yWallVelocity","Y-component of a velocity vector of the wall.");
  options.addConfigOption< CFreal >("xWallVelocity","X-component of a velocity vector of the wall.");
  options.addConfigOption< CFreal >("KWall","Wall value for turbulent intensity");
  options.addConfigOption< CFreal,Config::DynamicOption<> >("T","wall static temperature");
  options.addConfigOption< CFreal >("q","wall heat flux");
  options.addConfigOption< bool >("HeatFlux","bool to tell if the wall has constant heat flux (possibly zero), default true.");
  options.addConfigOption< CFuint,Config::DynamicOption<> >("ChangeToIsoT","Iteration after which to switch to an isothermal BC.");
}

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallTurb2D::BCNoSlipWallTurb2D(const std::string& name) :
  BCStateComputer(name),
  m_varSetTurb(CFNULL),
  m_diffVarTurb(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
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
  
  m_wallK = 1.e-8;
   setParameter("KWall",&m_wallK);
}

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallTurb2D::~BCNoSlipWallTurb2D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallTurb2D::computeGhostStates(const vector< State* >& intStates,
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
      m_ghostSolPhysData[EulerTerm::VX]  = 2.0*m_xWallVelocity-m_intSolPhysData[EulerTerm::VX];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::VY]  = 2.0*m_yWallVelocity-m_intSolPhysData[EulerTerm::VY];// negate velocity )-- average = 0
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
      m_ghostSolPhysData[EulerTerm::VX]  = 2.0*m_xWallVelocity-m_intSolPhysData[EulerTerm::VX];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::VY]  = 2.0*m_yWallVelocity-m_intSolPhysData[EulerTerm::VY];// negate velocity )-- average = 0
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
    
    if(nbTurbVars == 2)
    {
      //Compute distance to innerstate
      CFreal y0 = 1.e-2;//1.e-9;
    
      //avoid too small distances
      //y0 = std::max(y0, 10.e-10);
    
      const CFreal pdim =  m_ghostSolPhysData[EulerTerm::P] * m_varSetTurb->getModel()->getPressRef();
      const CFreal Tdim =  m_ghostSolPhysData[EulerTerm::T] * m_varSetTurb->getModel()->getTempRef();
      const CFreal mu = m_diffVarTurb->getModel().getDynViscosityDim(pdim, Tdim)/(m_diffVarTurb->getModel().getReferencePhysicalData())[NSTurbTerm::MU];
    
      CFreal nu = mu / m_ghostSolPhysData[EulerTerm::RHO];
    
      //this is not the best, but it avoids having to code another BC! because I
      //would have to dynamic cast to the KOmega varset to get the beta1
      const CFreal beta1 = 0.075;
      
      ///@todo here should this be adimensionalized (by the distance)???
      //Menter's definition
      const CFreal omegaWall = (10. * 6. * nu) / (beta1 * y0 * y0);
      
      cf_assert(omegaWall>0.0);
    
      const CFuint omegaID = iK+1; 

      m_ghostSolPhysData[omegaID] = omegaWall;
    }

    // set the ghost state from its physical data
    m_varSetTurb->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
    
    //CFLog(INFO, "ghostState: " << ghostState << "\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallTurb2D::computeGhostGradients
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
      tempGradG = tempGradI - 2.0*nTempGrad*normal + m_wallQ*normal;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallTurb2D::setup()
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

