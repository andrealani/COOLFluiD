#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/BCSubInletEulerVT2D.hh"

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
    BCSubInletEulerVT2D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNavierStokesModule >
  BCSubInletEulerVT2DProvider("SubInletEulerVT2D");

//////////////////////////////////////////////////////////////////////////////

void BCSubInletEulerVT2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Vx","x-velocity");
  options.addConfigOption< CFreal >("Vy","y-velocity");
  options.addConfigOption< CFreal >("T","temperature");
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletEulerVT2D::BCSubInletEulerVT2D(const std::string& name) :
  BCStateComputer(name),
  m_eulerVarSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_u(),
  m_v(),
  m_T()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_u = 0.0;
   setParameter("Vx",&m_u);

  m_v = 0.0;
   setParameter("Vy",&m_v);

  m_T = 0.0;
   setParameter("T",&m_T);
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletEulerVT2D::~BCSubInletEulerVT2D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletEulerVT2D::computeGhostStates(const vector< State* >& intStates,
                                                    vector< State* >& ghostStates,
                                                    const std::vector< RealVector >& normals,
                                                    const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  // get some data from the physical model
  const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
  const CFreal idGassConst = m_eulerVarSet->getModel()->getR();

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);

    cf_assert(intState.size() == 4);
    cf_assert(ghostState.size() == 4);

    // set the physical data starting from the inner state
    m_eulerVarSet->computePhysicalData(intState,m_intSolPhysData);

    //set the physical data for the ghost state
    m_ghostSolPhysData[EulerTerm::P]   = m_intSolPhysData[EulerTerm::P];
    m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::P]/(idGassConst*m_T);
    m_ghostSolPhysData[EulerTerm::VX]  = 2.0*m_u - m_intSolPhysData[EulerTerm::VX];
    m_ghostSolPhysData[EulerTerm::VY]  = 2.0*m_v - m_intSolPhysData[EulerTerm::VY];
    m_ghostSolPhysData[EulerTerm::V] = std::sqrt(m_ghostSolPhysData[EulerTerm::VX]*
				       m_ghostSolPhysData[EulerTerm::VX] +
				       m_ghostSolPhysData[EulerTerm::VY]*
				       m_ghostSolPhysData[EulerTerm::VY]);
    m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*m_ghostSolPhysData[EulerTerm::P]
				         + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
				         m_ghostSolPhysData[EulerTerm::V]*
				         m_ghostSolPhysData[EulerTerm::V])/m_ghostSolPhysData[EulerTerm::RHO];
    m_ghostSolPhysData[EulerTerm::A] = sqrt(gamma*m_ghostSolPhysData[EulerTerm::P]/
				       m_ghostSolPhysData[EulerTerm::RHO]);
    m_ghostSolPhysData[EulerTerm::T] = 2.0*m_T - m_intSolPhysData[EulerTerm::T];

    // set the ghost state from its physical data
    m_eulerVarSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletEulerVT2D::computeGhostGradients
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

void BCSubInletEulerVT2D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in BCSubInletEulerVT2D!");
  }

  // resize the physical data for internal and ghost solution points
  m_eulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_eulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );

  // non-dimensionalize pressure and temperature
  m_u /= m_eulerVarSet->getModel()->getVelRef ();
  m_v /= m_eulerVarSet->getModel()->getVelRef();
  m_T /= m_eulerVarSet->getModel()->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

