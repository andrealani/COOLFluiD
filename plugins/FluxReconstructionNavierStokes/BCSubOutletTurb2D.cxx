#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/BCSubOutletTurb2D.hh"

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
    BCSubOutletTurb2D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNavierStokesModule >
  BCSubOutletTurb2DProvider("BCSubOutletTurb2D");

//////////////////////////////////////////////////////////////////////////////

void BCSubOutletTurb2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("P","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

BCSubOutletTurb2D::BCSubOutletTurb2D(const std::string& name) :
  BCStateComputer(name),
  m_varSetTurb(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);
   
  m_pressure = 1.0;
  setParameter("P",&m_pressure);
}

//////////////////////////////////////////////////////////////////////////////

BCSubOutletTurb2D::~BCSubOutletTurb2D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubOutletTurb2D::computeGhostStates(const vector< State* >& intStates,
                                                    vector< State* >& ghostStates,
                                                    const std::vector< RealVector >& normals,
                                                    const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

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
    
    const CFreal pInner = m_intSolPhysData[EulerTerm::P];
    
    const CFreal pGhost = max(2.0*m_pressure - pInner,1e-10);

    //set the physical data for the ghost state
    m_ghostSolPhysData[EulerTerm::P]   = pGhost;
    m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::RHO];
    m_ghostSolPhysData[EulerTerm::VX]  = m_intSolPhysData[EulerTerm::VX];
    m_ghostSolPhysData[EulerTerm::VY]  = m_intSolPhysData[EulerTerm::VY];
    m_ghostSolPhysData[EulerTerm::V] = m_intSolPhysData[EulerTerm::V];
    m_ghostSolPhysData[EulerTerm::H] = (gammaDivGammaMinus1*m_ghostSolPhysData[EulerTerm::P]
				       + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
				       m_ghostSolPhysData[EulerTerm::V]*
				       m_ghostSolPhysData[EulerTerm::V])/m_ghostSolPhysData[EulerTerm::RHO];
    m_ghostSolPhysData[EulerTerm::A] = sqrt(gamma*m_ghostSolPhysData[EulerTerm::P]/
				       m_ghostSolPhysData[EulerTerm::RHO]);
    m_ghostSolPhysData[EulerTerm::T] = m_ghostSolPhysData[EulerTerm::P]/(R*m_ghostSolPhysData[EulerTerm::RHO]);
    m_ghostSolPhysData[EulerTerm::E] = m_ghostSolPhysData[EulerTerm::H] -
                                       (m_ghostSolPhysData[EulerTerm::P]/m_ghostSolPhysData[EulerTerm::RHO]);
    
    const CFuint iK = m_varSetTurb->getModel()->getFirstScalarVar(0);
    const CFuint nbTurbVars = m_varSetTurb->getModel()->getNbScalarVars(0);
    for(CFuint iTurb = 0; iTurb < nbTurbVars; iTurb++)
    {
      m_ghostSolPhysData[iK + iTurb] = m_intSolPhysData[iK + iTurb];
    }

    // set the ghost state from its physical data
    m_varSetTurb->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSubOutletTurb2D::computeGhostGradients
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
}

//////////////////////////////////////////////////////////////////////////////

void BCSubOutletTurb2D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;
  
  m_varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb2DVarSet>();
  m_varSetTurb->getModel()->resizePhysicalData(m_intSolPhysData);
  m_varSetTurb->getModel()->resizePhysicalData(m_ghostSolPhysData);

  m_pressure /= m_varSetTurb->getModel()->getPressRef();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

