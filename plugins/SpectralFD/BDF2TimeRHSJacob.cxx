#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSIdxMapping.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/BDF2TimeRHSJacob.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BDF2TimeRHSJacob,
                       SpectralFDMethodData,
                       SpectralFDModule>
BDF2_RHSJacob("BDF2_RHSJacob");

//////////////////////////////////////////////////////////////////////////////

BDF2TimeRHSJacob::BDF2TimeRHSJacob(const std::string& name) :
    PseudoSteadyStdTimeRHSJacob(name),
    socket_pastTimeRhs("pastTimeRhs")
{
}

//////////////////////////////////////////////////////////////////////////////

BDF2TimeRHSJacob::~BDF2TimeRHSJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void BDF2TimeRHSJacob::setup()
{
  PseudoSteadyStdTimeRHSJacob::setup();
}

//////////////////////////////////////////////////////////////////////////////

void BDF2TimeRHSJacob::addTimeResidual()
{
  // get factor for the residual (in the Jacobian)
  const CFreal resFactor = getMethodData().getResFactor();
  const CFreal resFactorP1 = resFactor + 1.;

  // get rhs and past states
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle<CFreal> pastTimeRhs = socket_pastTimeRhs.getDataHandle();

  // local (current processor) to global mapping
  const LSSIdxMapping& idxMapping =
      getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();

  // add time residual contribution (Backward Euler)
  const CFuint nbrSolPnts = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // get state
    const State& currState = *(*m_cellStates)[iSol];

    // get state ID
    const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

    // get past state
    const State& pastState = *pastStates[stateID];

    // add contribution to rhs and jacobian
    CFuint globalID = idxMapping.getColID(stateID)*m_nbrEqs;
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++globalID)
    {
      const CFreal dUdt = (currState[iEq] - pastState[iEq])*m_diagValues[iSol];
      rhs(stateID,iEq,m_nbrEqs) -= dUdt*resFactorP1 - resFactor*pastTimeRhs(stateID,iEq,m_nbrEqs);

      if (getMethodData().doComputeJacobian())
      {
        m_jacobMatrix->addValue(globalID, globalID, m_diagValues[iSol]*resFactorP1);
      }

      // store dUdt for next time step
      if(SubSystemStatusStack::getActive()->isLastStep())
      {
        pastTimeRhs(stateID,iEq,m_nbrEqs) = dUdt;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
BDF2TimeRHSJacob::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> >
      result = PseudoSteadyStdTimeRHSJacob::needsSockets();

  result.push_back(&socket_pastTimeRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
