#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/BCNoSlipWallHeatFluxNS3D.hh"
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
    BCNoSlipWallHeatFluxNS3D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNavierStokesModule >
  BCNoSlipWallHeatFluxNS3DProvider("NoSlipWallHeatFluxNS3D");

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallHeatFluxNS3D::BCNoSlipWallHeatFluxNS3D(const std::string& name) :
  BCStateComputer(name),
  m_eulerVarSet(CFNULL),
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

  m_legacyGhost = false;
   setParameter("LegacyGhost",&m_legacyGhost);
}

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallHeatFluxNS3D::~BCNoSlipWallHeatFluxNS3D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxNS3D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal,Config::DynamicOption<> >("T","wall static temperature");
  options.addConfigOption< CFreal,Config::DynamicOption<> >("q","wall heat flux");
  options.addConfigOption< bool >("HeatFlux","bool to tell if the wall has constant heat flux, default true.");
  options.addConfigOption< CFuint,Config::DynamicOption<> >("ChangeToIsoT","Iteration after which to switch to an isothermal BC.");
  options.addConfigOption< bool >("LegacyGhost","Use the previous ghost state (the wall state on an isothermal wall) "
    "instead of the reflected one (interior density, reversed velocity, reflected temperature), default false.");
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxNS3D::computeGhostStates(const vector< State* >& intStates,
                                                  vector< State* >& ghostStates,
                                                  const std::vector< RealVector >& normals,
                                                  const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());
  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  if (iter >= m_changeToIsoT && m_heatFlux)
  {
    m_heatFlux = false;
  }

  // get some physical data from the model
  const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intState   = (*intStates  [iState]);
    State& ghostState = (*ghostStates[iState]);

    cf_assert(intState.size() == 5);
    cf_assert(ghostState.size() == 5);

    // set the physical data starting from the inner state
    m_eulerVarSet->computePhysicalData(intState,m_intSolPhysData);
    
    if (!m_legacyGhost)
    {
      // interior density, velocity reversed (fixed wall), temperature reflected about the wall
      // temperature (isothermal wall) or copied (heat-flux wall), pressure from density and temperature
      const CFreal R = m_eulerVarSet->getModel()->getR();
      const CFreal innerT = m_intSolPhysData[EulerTerm::P]/(R*m_intSolPhysData[EulerTerm::RHO]);
      const CFreal ghostT = m_heatFlux ? innerT : max(2.0*m_wallT - innerT,0.01*m_wallT);
      const CFreal ghostRho = m_intSolPhysData[EulerTerm::RHO];
      const CFreal ghostP   = ghostRho*R*ghostT;

      m_ghostSolPhysData[EulerTerm::RHO] = ghostRho;
      m_ghostSolPhysData[EulerTerm::VX]  = -m_intSolPhysData[EulerTerm::VX];
      m_ghostSolPhysData[EulerTerm::VY]  = -m_intSolPhysData[EulerTerm::VY];
      m_ghostSolPhysData[EulerTerm::VZ]  = -m_intSolPhysData[EulerTerm::VZ];
      m_ghostSolPhysData[EulerTerm::V]   = m_intSolPhysData[EulerTerm::V];
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
      m_ghostSolPhysData[EulerTerm::VX]  = -m_intSolPhysData[EulerTerm::VX];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::VY]  = -m_intSolPhysData[EulerTerm::VY];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::VZ]  = -m_intSolPhysData[EulerTerm::VZ];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::V] = sqrt(m_ghostSolPhysData[EulerTerm::VX]*m_ghostSolPhysData[EulerTerm::VX]+m_ghostSolPhysData[EulerTerm::VY]*m_ghostSolPhysData[EulerTerm::VY]+m_ghostSolPhysData[EulerTerm::VZ]*m_ghostSolPhysData[EulerTerm::VZ]);
      m_ghostSolPhysData[EulerTerm::P]   = m_intSolPhysData[EulerTerm::P];
      m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*m_ghostSolPhysData[EulerTerm::P]
                                            + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
                                                  m_intSolPhysData[EulerTerm::V]*
                                                  m_intSolPhysData[EulerTerm::V]
                                         )/m_ghostSolPhysData[EulerTerm::RHO];
      m_ghostSolPhysData[EulerTerm::A] = sqrt(gamma*m_intSolPhysData[EulerTerm::P]/m_intSolPhysData[EulerTerm::RHO]);
      m_ghostSolPhysData[EulerTerm::T] = m_intSolPhysData[EulerTerm::T];
      m_ghostSolPhysData[EulerTerm::E] = m_ghostSolPhysData[EulerTerm::H] -
                                         (m_ghostSolPhysData[EulerTerm::P]/m_ghostSolPhysData[EulerTerm::RHO]);      
    }
    else
    {
      const CFreal R = m_eulerVarSet->getModel()->getR();
      const CFreal innerT = m_intSolPhysData[EulerTerm::P]/(R*m_intSolPhysData[EulerTerm::RHO]);
      const CFreal ghostT = m_wallT;
//      CFreal ghostP;
//      if (getMethodData().getUpdateVarStr() == "Cons")
//      {
//	ghostP = m_intSolPhysData[EulerTerm::RHO]*R*ghostT;
//      }
//      else
//      {
//	ghostP = m_intSolPhysData[EulerTerm::P];
//      }

      // set the physical data for the ghost state
      m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::P]/(R*ghostT);
      m_ghostSolPhysData[EulerTerm::VX]  = 0.; //-m_intSolPhysData[EulerTerm::VX];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::VY]  = 0.; //-m_intSolPhysData[EulerTerm::VY];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::VZ]  = 0.; //-m_intSolPhysData[EulerTerm::VZ];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::V] = sqrt(m_ghostSolPhysData[EulerTerm::VX]*m_ghostSolPhysData[EulerTerm::VX]+m_ghostSolPhysData[EulerTerm::VY]*m_ghostSolPhysData[EulerTerm::VY]+m_ghostSolPhysData[EulerTerm::VZ]*m_ghostSolPhysData[EulerTerm::VZ]);
      m_ghostSolPhysData[EulerTerm::P]   = m_intSolPhysData[EulerTerm::P]; //ghostP;
      m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*m_ghostSolPhysData[EulerTerm::P]
                                            + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
                                                  m_intSolPhysData[EulerTerm::V]*
                                                  m_intSolPhysData[EulerTerm::V]
                                         )/m_ghostSolPhysData[EulerTerm::RHO];
      m_ghostSolPhysData[EulerTerm::A] = sqrt(gamma*m_ghostSolPhysData[EulerTerm::P]/m_ghostSolPhysData[EulerTerm::RHO]);
      m_ghostSolPhysData[EulerTerm::T] = ghostT;
      m_ghostSolPhysData[EulerTerm::E] = m_ghostSolPhysData[EulerTerm::H] -
                                         (m_ghostSolPhysData[EulerTerm::P]/m_ghostSolPhysData[EulerTerm::RHO]);      
    }

    // set the ghost state from its physical data
    m_eulerVarSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxNS3D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                                     std::vector< std::vector< RealVector* > >& ghostGrads,
                                                     const std::vector< RealVector >& normals,
                                                     const std::vector< RealVector >& coords)
{
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    cf_assert(intGrads[iState].size() == 5);

    // normal
    const RealVector& normal = normals[iState];

    // pressure
//     RealVector& pGradI = *intGrads  [iState][0];
//     RealVector& pGradG = *ghostGrads[iState][0];
//     const CFreal nPGrad = pGradI[XX]*normal[XX] + pGradI[YY]*normal[YY] + pGradI[ZZ]*normal[ZZ];
//     pGradG = pGradI - 2.0*nPGrad*normal;
    *ghostGrads[iState][0] = *intGrads[iState][0];

    // velocity
    *ghostGrads[iState][1] = *intGrads[iState][1];
    *ghostGrads[iState][2] = *intGrads[iState][2];
    *ghostGrads[iState][3] = *intGrads[iState][3];
    
    if (m_heatFlux)
    {
      // temperature
      RealVector& tempGradI = *intGrads  [iState][4];
      RealVector& tempGradG = *ghostGrads[iState][4];
      const CFreal nTempGrad = tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY] + tempGradI[ZZ]*normal[ZZ] ;
      tempGradG = tempGradI - 2.0*nTempGrad*normal + m_wallQ*normal;
    }
    else
    {
      *ghostGrads[iState][4] = *intGrads[iState][4];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxNS3D::computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                                                  const std::vector< Framework::State* >& intStates,
                                                  const std::vector< Framework::State* >& ghostStates,
                                                  const std::vector< RealVector >& unitNormals,
                                                  const std::vector< RealVector >& flxPntCoords,
                                                  std::vector< RealVector* >& bndGradVars)
{
  const CFuint nbrStates = intStates.size();

  // g_b = a with zero velocity and, for an isothermal wall, the wall temperature
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    *bndGradVars[iState] = *gradVarsFlxPnt[iState];
    (*bndGradVars[iState])[1] = 0.;
    (*bndGradVars[iState])[2] = 0.;
    (*bndGradVars[iState])[3] = 0.;

    if (!m_heatFlux)
    {
      (*bndGradVars[iState])[4] = m_wallT;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxNS3D::computeBndStates(const std::vector< Framework::State* >& intStates,
                                                const std::vector< Framework::State* >& ghostStates,
                                                const std::vector< RealVector >& unitNormals,
                                                const std::vector< RealVector >& flxPntCoords,
                                                std::vector< RealVector* >& bndStates)
{
  const CFuint nbrStates = intStates.size();

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    m_eulerVarSet->computePhysicalData(*intStates[iState],m_intSolPhysData);

    // interior pressure, zero velocity, wall or interior temperature
    m_bndPrimState[0] = m_intSolPhysData[EulerTerm::P];
    m_bndPrimState[4] = m_heatFlux ? m_intSolPhysData[EulerTerm::T] : m_wallT;
    m_bndPrimState[1] = 0.;
    m_bndPrimState[2] = 0.;
    m_bndPrimState[3] = 0.;

    computeNSBoundaryState(*m_eulerVarSet,*intStates[iState],m_bndPrimState,*bndStates[iState]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxNS3D::computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                                               std::vector< std::vector< RealVector* > >& bndGrads,
                                               const std::vector< RealVector* >& bndStates,
                                               const std::vector< RealVector >& unitNormals,
                                               const std::vector< RealVector >& flxPntCoords)
{
  // q_b = q
  copyGradients(intGrads,bndGrads);
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxNS3D::constrainBndGrads(const RealVector& bndState,
                                                 std::vector< RealVector* >& bndGrads,
                                                 const RealVector& unitNormal,
                                                 const RealVector& flxPntCoord)
{
  if (m_heatFlux)
  {
    prescribeNSWallHeatFlux(*m_diffusiveVarSet,bndState,bndGrads,unitNormal,4,m_wallQ);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxNS3D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // get Euler 3D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSet in BCNoSlipWallHeatFluxNS3D!");
  }

  // resize the physical data for internal and ghost solution points
  m_eulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_eulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );

  // boundary primitive variables
  m_bndPrimState.resize(5);

  // diffusive variable set, for the prescribed heat flux
  m_diffusiveVarSet = getMethodData().getDiffusiveVar().d_castTo< NavierStokesVarSet >();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

