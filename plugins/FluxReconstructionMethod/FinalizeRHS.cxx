#include "Framework/MethodCommandProvider.hh"

#include "MathTools/MatrixInverter.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FinalizeRHS.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodCommandProvider<
    FinalizeRHS,FluxReconstructionSolverData,FluxReconstructionModule >
  FinalizeRHSProvider("StdFinalize");

//////////////////////////////////////////////////////////////////////////////

void FinalizeRHS::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

FinalizeRHS::FinalizeRHS(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_states("states"),
  socket_rhs("rhs"),
  m_updateToSolutionVecTrans(CFNULL),
  m_numJacob(CFNULL),
  m_inverter(CFNULL), 
  m_jacobDummy(),
  m_invJacobDummy(),
  m_state()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

}

//////////////////////////////////////////////////////////////////////////////

FinalizeRHS::~FinalizeRHS()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  FinalizeRHS::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FinalizeRHS::execute()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  RealVector tempRes(nbEqs);
  RealVector res(nbEqs);

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  
  for(CFuint iState = 0; iState < nbStates; ++iState) {

    // set and get the transformation matrix in the update variables
    const RealMatrix& matrix = computeNumericalTransMatrix(*states[iState]);

    // copy the rhs in a given temporary array
    const CFuint startID = iState*nbEqs;
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      tempRes[iEq] = rhs[startID + iEq];
    }

    // compute the transformed residual
    res = matrix*tempRes;
    
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      rhs[startID + iEq] = res[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& FinalizeRHS::computeNumericalTransMatrix(State& state)
{
  // this first transformed state HAS TO BE stored,
  // since the returned pointer will change pointee after
  // the second call to transform()
  m_state = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(&state));

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    // perturb the given component of the state vector
    m_numJacob->perturb(iVar, state[iVar]);

    const RealVector& tPertState = static_cast<RealVector&>
      (*m_updateToSolutionVecTrans->transform(&state));

    // compute the finite difference derivative of the flux
    m_numJacob->computeDerivative(m_state, tPertState, m_dState);

    m_jacobDummy.setColumn(m_dState,iVar);

    // restore the unperturbed value
    m_numJacob->restore(state[iVar]);
  }

  m_inverter->invert(m_jacobDummy, m_invJacobDummy);

  return m_invJacobDummy;
}


//////////////////////////////////////////////////////////////////////////////

void FinalizeRHS::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  FluxReconstructionSolverCom::setup();
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  
  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  m_state.resize(nbEqs,0.);
  m_dState.resize(nbEqs,0.);
  m_jacobDummy.resize(nbEqs, nbEqs, 0.);
  m_invJacobDummy.resize(nbEqs, nbEqs, 0.);
  
  m_inverter.reset(MatrixInverter::create(nbEqs, false));
  
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

