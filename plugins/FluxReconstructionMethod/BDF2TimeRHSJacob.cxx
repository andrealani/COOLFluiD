#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSIdxMapping.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/BlockAccumulator.hh"


#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/BDF2TimeRHSJacob.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BDF2TimeRHSJacob,
                       FluxReconstructionSolverData,
                       FluxReconstructionModule>
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
    // now the contribution to the jacobian matrix is calculated
    m_acc->setRowColIndex(iSol, (*m_cellStates)[iSol]->getLocalID());
  }
  
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // get state
    State *const currState = (*m_cellStates)[iSol];
    
   m_tempState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(currState));

    // get state ID
    const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

    // get past state
    const State *const pastState = m_updateToSolutionVecTrans->transform(pastStates[stateID]); //*pastStates[stateID];

    // add contribution to rhs and jacobian
    CFuint globalID = idxMapping.getColID(stateID)*m_nbrEqs;
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++globalID)
    {
      const CFreal dUdt = (m_tempState[iEq] - (*pastState)[iEq])*m_diagValues[iSol];
      rhs(stateID,iEq,m_nbrEqs) -= dUdt*resFactorP1 - resFactor*pastTimeRhs(stateID,iEq,m_nbrEqs);
      
      if((!getMethodData().isSysMatrixFrozen()) && getMethodData().doComputeJacobian())
      {
        // perturb the given component of the state vector
        m_numericalJacob->perturb(iEq, (*currState)[iEq]);

        const RealVector& tempPertState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(currState));

        // compute the finite difference derivative of the flux
        m_numericalJacob->computeDerivative(m_tempState,tempPertState,m_fluxDiff);

        // _fluxDiff corresponds to a column vector of the dU/dP matrix
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) 
        {
	  m_fluxDiff[iEq] *= m_diagValues[iSol]*resFactorP1;
        }

        m_acc->addValues(iSol, iSol, iEq, &m_fluxDiff[0]);

        // restore the unperturbed value
        m_numericalJacob->restore((*currState)[iEq]);

      }

//      if (getMethodData().doComputeJacobian())
//      {
//        m_jacobMatrix->addValue(globalID, globalID, m_diagValues[iSol]*resFactorP1);
//      }

      // store dUdt for next time step
      if((SubSystemStatusStack::getActive()->isLastStep())
       && (SubSystemStatusStack::getActive()->isSubIterationLastStep()))
      {
        pastTimeRhs(stateID,iEq,m_nbrEqs) = dUdt;
      }
    }
    
    if((!getMethodData().isSysMatrixFrozen()) && getMethodData().doComputeJacobian())
    {
      // add the values in the jacobian matrix
      //getMethodData().getLSSMatrix(0)->addValues(*m_acc);
      m_lss->getMatrix()->addValues(*m_acc);

      // reset to zero the entries in the block accumulator
      m_acc->reset();
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

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
