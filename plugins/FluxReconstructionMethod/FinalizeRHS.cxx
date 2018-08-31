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

FinalizeRHS::FinalizeRHS(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_states("states"),
  socket_rhs("rhs"),
  m_solPhysData(),
  m_jacobDummy(),
  m_invJacobDummy(),
  m_updateToSolutionVecTrans(CFNULL),
  m_solutionToUpdateMatTrans(CFNULL),
  m_state(),
  m_dState(),
  m_inverter(CFNULL),
  m_numJacob(CFNULL),
  m_tempRes()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);
  
  m_useAnalyticalMatrix = false;
  setParameter("useAnalyticalMatrix",&m_useAnalyticalMatrix);

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

void FinalizeRHS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >
    ("useAnalyticalMatrix", "Flag telling if to use analytical matrix."); 
}

//////////////////////////////////////////////////////////////////////////////

void FinalizeRHS::execute()
{
  // get the number of equations
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // get the state socket
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // get total nb of states and initialize vars
  const CFuint nbStates = states.size();

  // get rhs socket
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  
  // loop over states to transform residual
  for(CFuint iState = 0; iState < nbStates; ++iState) 
  {
    // set and get the transformation matrix in the update variables
    const RealMatrix& matrix = (m_useAnalyticalMatrix) ? computeAnalyticalTransMatrix(*states[iState]) : computeNumericalTransMatrix(*states[iState]);

    // copy the rhs in a given temporary array
    const CFuint startID = iState*nbEqs;
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) 
    {
      m_tempRes[iEq] = rhs[startID + iEq];
      if (abs(rhs[startID + iEq]) > 1.0e-3) CFLog(INFO,"eq: " << iEq << ", temp rhs: " << rhs[startID + iEq] << "\n");
    }

    // compute the transformed residual
    m_tempRes = matrix*m_tempRes;
    
    // store transformed residual in rhs
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) 
    {
      rhs[startID + iEq] = m_tempRes[iEq];
      if (abs(rhs[startID + iEq]) > 1.0e-3) CFLog(INFO,"eq: " << iEq << ", rhs: " << rhs[startID + iEq] << "\n");//states[iState]->getLocalID()==704
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
  if (state.getLocalID()==704) CFLog(VERBOSE, "origstate: " << state << ", transstate: " << m_state << ", coords: " << (state.getCoordinates()) << "\n");

  // get nb of equations
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // loop over vars
  for (CFuint iVar = 0; iVar < nbEqs; ++iVar) 
  {
    // perturb the given component of the state vector
    m_numJacob->perturb(iVar, state[iVar]);

    const RealVector& tPertState = static_cast<RealVector&> (*m_updateToSolutionVecTrans->transform(&state));

    // compute the finite difference derivative of the flux
    m_numJacob->computeDerivative(m_state, tPertState, m_dState);

    m_jacobDummy.setColumn(m_dState,iVar);

    // restore the unperturbed value
    m_numJacob->restore(state[iVar]);
  }

  // invert the Jacobian
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
  m_solutionToUpdateMatTrans = getMethodData().getSolToUpdateInUpdateMatTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  m_solutionToUpdateMatTrans->setup(2);
  
  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  m_state.resize(nbEqs,0.);
  m_dState.resize(nbEqs,0.);
  m_jacobDummy.resize(nbEqs, nbEqs, 0.);
  m_invJacobDummy.resize(nbEqs, nbEqs, 0.);
  
  m_inverter.reset(MatrixInverter::create(nbEqs, false));
  
  m_tempRes.resize(nbEqs);
}

//////////////////////////////////////////////////////////////////////////////

void FinalizeRHS::unsetup()
{
  CFAUTOTRACE;

  // unsetup of the parent class
  FluxReconstructionSolverCom::unsetup();
  
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

