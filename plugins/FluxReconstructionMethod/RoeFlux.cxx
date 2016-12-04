#include <iterator>

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/CFSide.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/RoeFlux.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    RoeFlux,FluxReconstructionSolverData,RiemannFlux,FluxReconstructionModule >
  RoeFluxProvider("RoeFlux");

//////////////////////////////////////////////////////////////////////////////

RoeFlux::RoeFlux(const std::string& name) :
  RiemannFlux(name),
  m_linearizer(),
  m_solutionToLinearVarTrans(),
  m_updateToSolutionVarTrans(),
  m_sumFlux(),
  m_rightEv(),
  m_leftEv(),
  m_eValues(),
  m_rightEvalues(),
  m_leftEvalues(),
  m_absEvalues(),
  m_absJacob(),
  m_linStates()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

RoeFlux::~RoeFlux()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& RoeFlux::computeFlux(State& lState,
                                 State& rState,
                                 const RealVector& normal)
{
  SafePtr< ConvectiveVarSet > updateVarSet = getMethodData().getUpdateVar();

  // store the left and right states in vector
  m_updateStates[LEFT ] = &lState;
  m_updateStates[RIGHT] = &rState;

  // compute physical data for the left and the right internal flux points
  updateVarSet->computePhysicalData(lState, m_pData[LEFT]);
  updateVarSet->computePhysicalData(rState, m_pData[RIGHT]);
  
  // flux for right and left state (the physical data must be passed here!)
  m_sumFlux  = updateVarSet->getFlux()(m_pData[LEFT], normal);
  m_sumFlux += updateVarSet->getFlux()(m_pData[RIGHT], normal);
  
  // transform from update states (which are stored) to solution states (in which the equations are written)
  m_solStates[LEFT ] = getMethodData().getUpdateToSolutionVecTrans()->transform(m_updateStates[LEFT ]);              
  m_solStates[RIGHT] = getMethodData().getUpdateToSolutionVecTrans()->transform(m_updateStates[RIGHT]);  
  
  // transform from solution states to linear states
  m_linStates[LEFT ] = m_solutionToLinearVarTrans->transform(m_solStates[LEFT ]);              
  m_linStates[RIGHT] = m_solutionToLinearVarTrans->transform(m_solStates[RIGHT]);

  // linearize the states (for instance: compute the Roe averaged values) AND set the physical data of the model! This is done inside this function
  m_linearizer->linearize(m_linStates);

  // set the eigenvectors and eigenvalues of the linearized jacobian
  // USING THE PHYSICAL DATA THAT WAS SET IN linearize();
  updateVarSet->computeEigenValuesVectors(m_rightEv, m_leftEv, m_eValues, normal);

  // set the abs of the  eigen values (the implementation of this
  // function change if there are entropy or carbuncle fixes)
  SetAbsEigenValues();

  // abs of the jacobian
  m_absJacob = m_rightEv*(m_absEvalues*m_leftEv);

  // compute the Riemann flux
  m_rFlux = 0.5*(m_sumFlux - m_absJacob*(rState - lState));

  return m_rFlux;
}

//////////////////////////////////////////////////////////////////////////////

    
RealVector& RoeFlux::computeFlux(State& lState,RealVector& lExtraVars,
                                 State& rState,RealVector& rExtraVars,
                                 const RealVector& normal)
{
  // There is no implementation for extravars yet.
  return computeFlux(lState,rState,normal);
}

//////////////////////////////////////////////////////////////////////////////

    
void RoeFlux::SetAbsEigenValues()
{
  m_absEvalues = abs(m_eValues);
}

//////////////////////////////////////////////////////////////////////////////

void RoeFlux::setup()
{
  CFAUTOTRACE;

  RiemannFlux::setup();

  // resize variables
  m_rightEv.resize(m_nbrEqs,m_nbrEqs);
  m_leftEv.resize(m_nbrEqs,m_nbrEqs);
  m_eValues.resize(m_nbrEqs);
  m_rightEvalues.resize(m_nbrEqs);
  m_leftEvalues.resize(m_nbrEqs);
  m_absEvalues.resize(m_nbrEqs);
  m_absJacob.resize(m_nbrEqs,m_nbrEqs);
  m_sumFlux.resize(m_nbrEqs);

  // get the name of the physical model
  std::string physicsName = PhysicalModelStack::getActive()->getImplementor()->getConvectiveName();

  // get the varset names
  std::string updateVarName = getMethodData().getUpdateVarStr();
  std::string solutionVarName = getMethodData().getSolutionVarStr();
  std::string linearVarName = getMethodData().getLinearVarStr();

  // create the linearizer
  std::string linearizerName = physicsName + "Linear" + linearVarName;
  CFLog(INFO, "RoeFlux::setup() => linearizerName = " << linearizerName << "\n");
  m_linearizer = Environment::Factory<JacobianLinearizer>::getInstance().
    getProvider(linearizerName)->create(PhysicalModelStack::getActive());
  cf_assert(m_linearizer.isNotNull());

  // create the solution to linear variables transformer
  std::string solutionToLinearVarName
    = VarSetTransformer::getProviderName(physicsName,solutionVarName,linearVarName);
  CFLog(INFO, "RoeFlux::setup() => solutionToLinearVarName = " << solutionToLinearVarName << "\n");
  m_solutionToLinearVarTrans =
    Environment::Factory<VarSetTransformer>::getInstance().getProvider
    (solutionToLinearVarName)->create(PhysicalModelStack::getActive()->getImplementor());
  cf_assert(m_solutionToLinearVarTrans.isNotNull());

  // create the update to solution variable transformer
  std::string updateToSolutionVarName
    = VarSetTransformer::getProviderName(physicsName,updateVarName,solutionVarName);
  CFLog(INFO, "RoeFlux::setup() => updateToSolutionVarName = " << updateToSolutionVarName << "\n");
  m_updateToSolutionVarTrans =
    Environment::Factory<VarSetTransformer>::getInstance().getProvider
    (updateToSolutionVarName)->create(PhysicalModelStack::getActive()->getImplementor());
  cf_assert(m_updateToSolutionVarTrans.isNotNull());

  //
  m_linearizer->setMaxNbStates(2);
  m_solutionToLinearVarTrans->setup(2*m_maxNbrFlxPnts);
  m_updateToSolutionVarTrans->setup(2*m_maxNbrFlxPnts);
  m_linStates.resize(2);

}

//////////////////////////////////////////////////////////////////////////////

void RoeFlux::unsetup()
{
  CFAUTOTRACE;

  RiemannFlux::unsetup();

}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
