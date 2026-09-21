#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/BCNoSlipWallTurb3D.hh"
#include "FluxReconstructionNavierStokes/NSBoundaryState.hh"

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
    BCNoSlipWallTurb3D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNavierStokesModule >
  BCNoSlipWallTurb3DProvider("BCNoSlipWallTurb3D");

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallTurb3D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("yWallVelocity","Y-component of a velocity vector of the wall.");
  options.addConfigOption< CFreal >("xWallVelocity","X-component of a velocity vector of the wall.");
  options.addConfigOption< CFreal >("zWallVelocity","Z-component of a velocity vector of the wall.");
  options.addConfigOption< CFreal >("KWall","Wall value for turbulent intensity");
  options.addConfigOption< CFreal,Config::DynamicOption<> >("T","wall static temperature");
  options.addConfigOption< CFreal >("q","wall heat flux");
  options.addConfigOption< bool >("HeatFlux","bool to tell if the wall has constant heat flux (possibly zero), default true.");
  options.addConfigOption< CFuint,Config::DynamicOption<> >("ChangeToIsoT","Iteration after which to switch to an isothermal BC.");
  options.addConfigOption< bool >("LegacyGhost","Use the previous ghost state (the wall state) "
    "instead of the reflected one (interior density, reflected velocity and temperature), default false.");
  options.addConfigOption< CFreal >("WallDist","Characteristic distance of first sol pnt from the wall.");
  options.addConfigOption< CFreal >("OmegaWallFactor","Factor by which to multiply omegaWall each iteration until it is the theoretical value (Default 1.01).");
  options.addConfigOption< CFuint >("ImposeOmegaWallIter","Iteration at which to impose theoretical omegaWall value.");
}

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallTurb3D::BCNoSlipWallTurb3D(const std::string& name) :
  BCStateComputer(name),
  m_varSetTurb(CFNULL),
  m_diffVarTurb(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_prevIter(0)
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

  m_legacyGhost = false;
   setParameter("LegacyGhost",&m_legacyGhost);
   
  m_xWallVelocity = 0.0;
   setParameter("xWallVelocity",&m_xWallVelocity);
  
  m_yWallVelocity = 0.0;
   setParameter("yWallVelocity",&m_yWallVelocity);
   
  m_zWallVelocity = 0.0;
  setParameter("zWallVelocity",&m_zWallVelocity);
  
  m_wallK = 0.0;
   setParameter("KWall",&m_wallK);
   
   m_wallDist = 1.0e-5;
   setParameter("WallDist",&m_wallDist);
   
   m_imposeOmegaWallIter = 0;
   setParameter("ImposeOmegaWallIter",&m_imposeOmegaWallIter);
   
   m_omegaWallFactor = 1.01;
   setParameter("OmegaWallFactor",&m_omegaWallFactor);
}

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallTurb3D::~BCNoSlipWallTurb3D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallTurb3D::computeGhostStates(const vector< State* >& intStates,
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
    
    if (!m_legacyGhost)
    {
      // interior density, velocity reflected about the wall velocity, temperature reflected about the
      // wall temperature (isothermal wall) or copied (heat-flux wall), pressure from density and temperature
      const CFreal innerT = m_intSolPhysData[EulerTerm::P]/(R*m_intSolPhysData[EulerTerm::RHO]);
      const CFreal ghostT = m_heatFlux ? innerT : max(2.0*m_wallT - innerT,0.01*m_wallT);
      const CFreal ghostRho = m_intSolPhysData[EulerTerm::RHO];
      const CFreal ghostP   = ghostRho*R*ghostT;

      m_ghostSolPhysData[EulerTerm::RHO] = ghostRho;
      m_ghostSolPhysData[EulerTerm::VX]  = 2.0*m_xWallVelocity - m_intSolPhysData[EulerTerm::VX];
      m_ghostSolPhysData[EulerTerm::VY]  = 2.0*m_yWallVelocity - m_intSolPhysData[EulerTerm::VY];
      m_ghostSolPhysData[EulerTerm::VZ]  = 2.0*m_zWallVelocity - m_intSolPhysData[EulerTerm::VZ];
      m_ghostSolPhysData[EulerTerm::V]   = sqrt(m_ghostSolPhysData[EulerTerm::VX]*m_ghostSolPhysData[EulerTerm::VX] +
                                              m_ghostSolPhysData[EulerTerm::VY]*m_ghostSolPhysData[EulerTerm::VY] +
                                              m_ghostSolPhysData[EulerTerm::VZ]*m_ghostSolPhysData[EulerTerm::VZ]);
      m_ghostSolPhysData[EulerTerm::P]   = ghostP;
      m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*ghostP + 0.5*ghostRho*
                                            m_ghostSolPhysData[EulerTerm::V]*m_ghostSolPhysData[EulerTerm::V])/ghostRho;
      m_ghostSolPhysData[EulerTerm::A]   = sqrt(gamma*ghostP/ghostRho);
      m_ghostSolPhysData[EulerTerm::T]   = ghostT;
      m_ghostSolPhysData[EulerTerm::E]   = m_ghostSolPhysData[EulerTerm::H] - ghostP/ghostRho;
    }
    else if (m_heatFlux)
    {
      // set the physical data for the ghost state
      m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::RHO];
      m_ghostSolPhysData[EulerTerm::VX]  = m_xWallVelocity;//2.0*m_xWallVelocity-m_intSolPhysData[EulerTerm::VX];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::VY]  = m_yWallVelocity;//2.0*m_yWallVelocity-m_intSolPhysData[EulerTerm::VY];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::VZ]  = m_zWallVelocity;//2.0*m_zWallVelocity-m_intSolPhysData[EulerTerm::VZ];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::V] = sqrt(m_ghostSolPhysData[EulerTerm::VX]*m_ghostSolPhysData[EulerTerm::VX]+m_ghostSolPhysData[EulerTerm::VY]*m_ghostSolPhysData[EulerTerm::VY]+m_ghostSolPhysData[EulerTerm::VZ]*m_ghostSolPhysData[EulerTerm::VZ]);
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
      m_ghostSolPhysData[EulerTerm::VZ]  = m_zWallVelocity;//2.0*m_zWallVelocity-m_intSolPhysData[EulerTerm::VZ];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::V] = sqrt(m_ghostSolPhysData[EulerTerm::VX]*m_ghostSolPhysData[EulerTerm::VX]+m_ghostSolPhysData[EulerTerm::VY]*m_ghostSolPhysData[EulerTerm::VY]+m_ghostSolPhysData[EulerTerm::VZ]*m_ghostSolPhysData[EulerTerm::VZ]);
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

    // wall values for the wall model of the turbulence variables: the ghost itself with the legacy
    // ghost, otherwise the interior pressure with the wall (isothermal) or interior (heat flux) temperature
    const CFreal wallP   = m_legacyGhost ? m_ghostSolPhysData[EulerTerm::P] : m_intSolPhysData[EulerTerm::P];
    const CFreal wallT   = m_legacyGhost ? m_ghostSolPhysData[EulerTerm::T] :
                           (m_heatFlux ? m_intSolPhysData[EulerTerm::T] : m_wallT);
    const CFreal wallRho = m_legacyGhost ? m_ghostSolPhysData[EulerTerm::RHO] : wallP/(R*wallT);

    m_ghostSolPhysData[iK] = m_wallK;
    
    // check if it is k-omega and not SA
    if(nbTurbVars == 2 || nbTurbVars == 4)
    {
      //Compute distance to innerstate
      CFreal y0 = m_wallDist;//1.e-9;
    
      //avoid too small distances
      //y0 = std::max(y0, 10.e-10);
    
      const CFreal pdim =  wallP * m_varSetTurb->getModel()->getPressRef();
      const CFreal Tdim =  wallT * m_varSetTurb->getModel()->getTempRef();
      const CFreal mu = m_diffVarTurb->getModel().getDynViscosityDim(pdim, Tdim)/(m_diffVarTurb->getModel().getReferencePhysicalData())[NSTurbTerm::MU];
    
      CFreal nu = mu / wallRho;
    
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
      m_ghostSolPhysData[EulerTerm::GAMMA] = m_intSolPhysData[EulerTerm::GAMMA];

      // gamma
      m_ghostSolPhysData[iK+2] = m_intSolPhysData[iK+2];
      
      // Re
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

void BCNoSlipWallTurb3D::computeGhostGradients
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
      RealVector& tempGradI = *intGrads  [iState][4];
      RealVector& tempGradG = *ghostGrads[iState][4];
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
  if (nbrGradVars > 7)
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

void BCNoSlipWallTurb3D::computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                                                  const std::vector< Framework::State* >& intStates,
                                                  const std::vector< Framework::State* >& ghostStates,
                                                  const std::vector< RealVector >& unitNormals,
                                                  const std::vector< RealVector >& flxPntCoords,
                                                  std::vector< RealVector* >& bndGradVars)
{
  const CFuint nbrStates = intStates.size();
  // index of k in a state and in the gradient variables, after p, the velocity and T; not the
  // index of k in the physical data, which is what getFirstScalarVar returns
  const CFuint iK = DIM_3D+2;
  const CFuint nbTurbVars = m_varSetTurb->getModel()->getNbScalarVars(0);

  // g_b = a with the wall values: the wall velocity, the wall temperature of an isothermal
  // wall, and KWall and the wall value of log-omega the ghost state carries
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const RealVector& ghostState = *ghostStates[iState];
    RealVector& bndGradVarsState = *bndGradVars[iState];

    bndGradVarsState = *gradVarsFlxPnt[iState];
    bndGradVarsState[1] = m_xWallVelocity;
    bndGradVarsState[2] = m_yWallVelocity;
    bndGradVarsState[3] = m_zWallVelocity;
    if (!m_heatFlux)
    {
      bndGradVarsState[DIM_3D+1] = m_wallT;
    }
    bndGradVarsState[iK] = ghostState[iK];
    if (nbTurbVars == 2 || nbTurbVars == 4)
    {
      bndGradVarsState[iK+1] = ghostState[iK+1];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallTurb3D::computeBndStates(const std::vector< Framework::State* >& intStates,
                                                const std::vector< Framework::State* >& ghostStates,
                                                const std::vector< RealVector >& unitNormals,
                                                const std::vector< RealVector >& flxPntCoords,
                                                std::vector< RealVector* >& bndStates)
{
  const CFuint nbrStates = intStates.size();

  // U_b is the wall state: the interior pressure, the wall velocity, the wall temperature of an
  // isothermal wall or the interior one of a heat-flux wall, and the wall values of the
  // turbulence variables the ghost state carries
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const RealVector& intState = *intStates[iState];
    RealVector& bndState = *bndStates[iState];

    bndState = *ghostStates[iState];
    bndState[0] = intState[0];
    bndState[1] = m_xWallVelocity;
    bndState[2] = m_yWallVelocity;
    bndState[3] = m_zWallVelocity;
    bndState[DIM_3D+1] = m_heatFlux ? intState[DIM_3D+1] : m_wallT;
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallTurb3D::computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                                               std::vector< std::vector< RealVector* > >& bndGrads,
                                               const std::vector< RealVector* >& bndStates,
                                               const std::vector< RealVector >& unitNormals,
                                               const std::vector< RealVector >& flxPntCoords)
{
  // q_b = q
  copyGradients(intGrads,bndGrads);

  // no diffusive flux of gamma and Re_theta through the wall: their normal component is removed
  const CFuint nbTurbVars = m_varSetTurb->getModel()->getNbScalarVars(0);
  if (nbTurbVars == 4)
  {
    // index of k in the gradients, as in computeBndGradVars
    const CFuint iK = DIM_3D+2;
    const CFuint nbrFlxPnts = intGrads.size();
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      removeNormalComponent(*bndGrads[iFlx][iK+2],unitNormals[iFlx]);
      removeNormalComponent(*bndGrads[iFlx][iK+3],unitNormals[iFlx]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallTurb3D::constrainBndGrads(const RealVector& bndState,
                                                 std::vector< RealVector* >& bndGrads,
                                                 const RealVector& unitNormal,
                                                 const RealVector& flxPntCoord)
{
  if (m_heatFlux)
  {
    prescribeNSWallHeatFlux(*m_diffVarTurb,bndState,bndGrads,unitNormal,DIM_3D+1,m_wallQ);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallTurb3D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  m_varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb3DVarSet>();

  m_diffVarTurb = getMethodData().getDiffusiveVar().d_castTo<DiffTurb3DVarSet>();

  m_varSetTurb->getModel()->resizePhysicalData(m_intSolPhysData);
  m_varSetTurb->getModel()->resizePhysicalData(m_ghostSolPhysData);
  
  m_xWallVelocity /= m_varSetTurb->getModel()->getVelRef();
  m_yWallVelocity /= m_varSetTurb->getModel()->getVelRef();
  m_zWallVelocity /= m_varSetTurb->getModel()->getVelRef();
  
  cf_assert(m_wallK >= 0.0);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

