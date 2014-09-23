#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/ComputeDiagBlockJacobMatrByPert.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeDiagBlockJacobMatrByPert, LUSGSIteratorData, LUSGSMethodModule>
    computeDiagBlockJacobMatrByPertProvider("DiagBlockJacobMatrByPert");

//////////////////////////////////////////////////////////////////////////////

ComputeDiagBlockJacobMatrByPert::ComputeDiagBlockJacobMatrByPert(std::string name) :
  LUSGSIteratorCom(name),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  socket_isStatesSetParUpdatable("isStatesSetParUpdatable"),
  socket_rhsCurrStatesSet("rhsCurrStatesSet"),
  socket_statesSetIdx("statesSetIdx"),
  socket_states("states"),
  socket_statesSetStateIDs("statesSetStateIDs"),
  m_numJacob(CFNULL),
  m_unpertRhs(CFNULL),
  m_nbrEqs(),
  m_stateIdx(),
  m_eqIdx()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiagBlockJacobMatrByPert::execute()
{
  CFAUTOTRACE;

  if (getMethodData().stopStatesLoop() && getMethodData().stopEqsLoop())
  // store unperturbed spacetime residual (rhs) for one states set
  {
    // Gets the rhs vectors datahandle
    DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();

    // Derefence the auxiliary rhs variable
    RealVector& unpertRhs = *m_unpertRhs;
    cf_assert(rhsCurrStatesSet.size() == unpertRhs.size());

    // copy the unperturbed spacetime residual (rhs)
    const CFuint nbrRes = rhsCurrStatesSet.size();
    for (CFuint iRes = 0; iRes < nbrRes; ++iRes)
    {
      unpertRhs[iRes] = rhsCurrStatesSet[iRes];
    }

    // clear states and equation index
    m_stateIdx = 0;
    m_eqIdx = 0;
  }
  else if (getMethodData().beforePertResComputation())
  {
    // Gets current states set index
    DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();
    const CFuint currStatesSetIdx = statesSetIdx[0];

    // set second element to 0 --> update coefficient is not recomputed
    statesSetIdx[1] = 0;

    // get isStatesSetParUpdatable data handle
    DataHandle< bool > isStatesSetParUpdatable = socket_isStatesSetParUpdatable.getDataHandle();

    if (isStatesSetParUpdatable[currStatesSetIdx])
    {
      // Gets the current states IDs
      DataHandle< vector< CFuint > > statesSetStateIDs = socket_statesSetStateIDs.getDataHandle();
      vector< CFuint > currStatesIDs = statesSetStateIDs[currStatesSetIdx];

      // Gets state datahandle
      DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

      // derefence the state we are currently perturbing
      cf_assert(m_stateIdx < currStatesIDs.size());
      const CFuint pertStateID = currStatesIDs[m_stateIdx];
      State& pertState = *states[pertStateID];

      // perturb physical variable in state
      m_numJacob->perturb(m_eqIdx,pertState[m_eqIdx]);
    }
  }
  else
  {
    // Gets current states set index
    DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();
    const CFuint currStatesSetIdx = statesSetIdx[0];

    // set second element to 1 --> update coefficient is recomputed
    statesSetIdx[1] = 1;

    // get isStatesSetParUpdatable data handle
    DataHandle< bool > isStatesSetParUpdatable = socket_isStatesSetParUpdatable.getDataHandle();

    // Gets the current states IDs
    DataHandle< vector< CFuint > > statesSetStateIDs = socket_statesSetStateIDs.getDataHandle();
    vector< CFuint > currStatesIDs = statesSetStateIDs[currStatesSetIdx];
    const CFuint nbrStates = currStatesIDs.size();

    if (isStatesSetParUpdatable[currStatesSetIdx])
    {
      // Gets state datahandle
      DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

      // derefence the state we are currently perturbing
      const CFuint pertStateID = currStatesIDs[m_stateIdx];
      State& pertState = *states[pertStateID];

      // restore physical variable in state
      m_numJacob->restore(pertState[m_eqIdx]);

      // Gets the current diagonal matrix
      DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();
      cf_assert(diagBlockJacobMatr.size() == getMethodData().getNbrStatesSets());
      RealMatrix& diagMatr = diagBlockJacobMatr[currStatesSetIdx];
      const CFuint resSize = diagMatr.nbRows();

      // Gets the rhs vectors datahandle
      DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();

      // Derefence the auxiliary rhs variable
      RealVector& unpertRhs = *m_unpertRhs;
      cf_assert(rhsCurrStatesSet.size() == unpertRhs.size());

      // vector for the perturbed rhs and the derivative of the rhs
      RealVector pertRhs (unpertRhs.size());
      RealVector derivRhs(unpertRhs.size());

      // copy perturbed rhs
      for (CFuint iRes = 0; iRes < resSize; ++iRes)
      {
        pertRhs[iRes] = rhsCurrStatesSet[iRes];
      }

      // compute the finite difference derivative
      m_numJacob->computeDerivative(pertRhs,unpertRhs,derivRhs);

      // add the derivatives to the current jacobian
      const CFuint colIdx = m_stateIdx*m_nbrEqs + m_eqIdx;
      for (CFuint iRes = 0; iRes < resSize; ++iRes)
      {
        diagMatr(iRes,colIdx) = derivRhs[iRes];
      }
    }

    // update state and equation indexes
    ++m_eqIdx;
    if (m_eqIdx == m_nbrEqs)
    {
      m_eqIdx = 0;
      getMethodData().setStopEqsLoop(true);
      ++m_stateIdx;
      if (m_stateIdx == nbrStates)
      {
        m_stateIdx = 0;
        getMethodData().setStopStatesLoop(true);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiagBlockJacobMatrByPert::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  LUSGSIteratorCom::setup();

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get the auxiliary rhs variable
  m_unpertRhs = getMethodData().getResAux();

  // get the number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > ComputeDiagBlockJacobMatrByPert::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_diagBlockJacobMatr     );
  result.push_back(&socket_isStatesSetParUpdatable);
  result.push_back(&socket_statesSetIdx);
  result.push_back(&socket_rhsCurrStatesSet);
  result.push_back(&socket_states);
  result.push_back(&socket_statesSetStateIDs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
