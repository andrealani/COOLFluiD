#include <iterator>

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/CFSide.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/RoeFlux.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    RoeFlux,SpectralFVMethodData,RiemannFlux,SpectralFVModule >
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
  m_updateStates(),
  m_updateStatesInFlxPnt(),
  m_extraVarsInFlxPnt(),
  m_solStates(),
  m_solStatesInFlxPnt(),
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

RealVector& RoeFlux::computeFlux(State& lState, RealVector& lExtraVars,
                                 State& rState, RealVector& rExtraVars,
                                 const RealVector& normal)
{
  // store the left and right states in vector
  m_updateStatesInFlxPnt[LEFT ] = &lState;
  m_updateStatesInFlxPnt[RIGHT] = &rState;
  m_extraVarsInFlxPnt   [LEFT ] = &lExtraVars;
  m_extraVarsInFlxPnt   [RIGHT] = &rExtraVars;

 // get the update VarSet
  SafePtr< ConvectiveVarSet > updateVarSet = getMethodData().getUpdateVar();

  // transform from update states (which are stored) to solution states (in which the equations are written)
 /// @todo broken after release 2009.3
//   m_solStates = m_updateToSolutionVarTrans->transformFromRef(&m_updateStatesInFlxPnt);

  // transform from solution states to linear states
  m_linStates = m_solutionToLinearVarTrans->transform(m_solStates);

  // linearize the states (for instance: compute the Roe averaged values) AND set the physical data of the model! This is done inside this function
  m_linearizer->linearize(*m_linStates);

  // set the eigenvectors and eigenvalues of the linearized jacobian
  // USING THE PHYSICAL DATA THAT WAS SET IN linearize();
  updateVarSet->computeEigenValuesVectors(m_rightEv, m_leftEv, m_eValues, normal);

  // set the abs of the  eigen values (the implementation of this
  // function change if there are entropy or carbuncle fixes)
  SetAbsEigenValues();

  // abs of the jacobian
  m_absJacob = m_rightEv*(m_absEvalues*m_leftEv);

  // compute the states data associated with lState and rState (must be called before computing the flux)
 /// @todo broken after release 2009.3
//   updateVarSet->computeStatesData(m_updateStatesInFlxPnt,m_extraVarsInFlxPnt,2);

  // flux for the right and left state (the UPDATE variables have to be passed here, this is implicitly assumed!!!
  m_sumFlux  = updateVarSet->getFlux()(lState, normal);
  m_sumFlux += updateVarSet->getFlux()(rState, normal);

  // compute the Riemann flux
  m_rFlux = 0.5*(m_sumFlux - m_absJacob*(rState - lState));

  return m_rFlux;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& RoeFlux::computeNumDamping(State& lState, RealVector& lExtraVars,
                                       State& rState, RealVector& rExtraVars,
                                       const RealVector& normal)
{
  // store the left and right states in vector
  m_updateStates[LEFT ] = &lState;
  m_updateStates[RIGHT] = &rState;

 // get the update VarSet
  SafePtr< ConvectiveVarSet > updateVarSet = getMethodData().getUpdateVar();

  // transform from update states (which are stored) to solution states (in which the equations are written)
 /// @todo broken after release 2009.3
//   m_solStates = m_updateToSolutionVarTrans->transformFromRef(&m_updateStates);

  // transform from solution states to linear states
  m_linStates = m_solutionToLinearVarTrans->transform(m_solStates);

  // linearize the states (for instance: compute the Roe averaged values) AND set the physical data of the model! This is done inside this function
  m_linearizer->linearize(*m_linStates);

  // set the eigenvectors and eigenvalues of the linearized jacobian
  // USING THE PHYSICAL DATA THAT WAS SET IN linearize();
  updateVarSet->computeEigenValuesVectors(m_rightEv, m_leftEv, m_eValues, normal);

  // set the abs of the  eigen values (the implementation of this
  // function change if there are entropy or carbuncle fixes)
  SetAbsEigenValues();

  // abs of the jacobian
  m_absJacob = m_rightEv*(m_absEvalues*m_leftEv);

  m_numDamping = m_absJacob*(rState - lState);

  return m_numDamping;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& RoeFlux::computeFlux(vector< Framework::State* >& lState,
                                           vector< Framework::State* >& rState,
                                           const vector< RealVector >& normal,
                                           const CFuint nbrFlxPnts)
{
  cf_assert(nbrFlxPnts <= lState.size());
  cf_assert(nbrFlxPnts <= rState.size());
  cf_assert(nbrFlxPnts <= normal.size());

  // store the left and right states in vector
  CFuint flxID = 0;
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_updateStates[flxID] = lState[iFlx]; ++flxID;
    m_updateStates[flxID] = rState[iFlx]; ++flxID;
  }

  // get the update VarSet
  SafePtr< ConvectiveVarSet > updateVarSet = getMethodData().getUpdateVar();

  // transform from update states (which are stored) to solution states (in which the equations are written)
 /// @todo broken after release 2009.3
//   m_solStates = m_updateToSolutionVarTrans->transformFromRef(&m_updateStates);

  // compute the Riemann flux for each flux point
  vector< Framework::State* >::iterator solStateItr = m_solStates->begin();
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    // last position to be copied
    vector< Framework::State* >::iterator lastSolState = solStateItr + 2;

    // copy the (pointers to) left and right states to auxiliary variable
    copy(solStateItr,lastSolState,m_solStatesInFlxPnt.begin());

    // transform from solution states to linear states
    m_linStates = m_solutionToLinearVarTrans->transform(&m_solStatesInFlxPnt);

    // linearize the states (for instance: compute the Roe averaged values) AND set the physical data of the model! This is done inside this function
    m_linearizer->linearize(*m_linStates);

    // set the eigenvectors and eigenvalues of the linearized jacobian
    // USING THE PHYSICAL DATA THAT WAS SET IN linearize();
    updateVarSet->computeEigenValuesVectors(m_rightEv, m_leftEv, m_eValues, normal[iFlx]);

    // set the abs of the  eigen values (the implementation of this
    // function change if there are entropy or carbuncle fixes)
    SetAbsEigenValues();

    // abs of the jacobian
    m_absJacob = m_rightEv*(m_absEvalues*m_leftEv);

    // dereference left and right update states
    State& lUpdState = *lState[iFlx];
    State& rUpdState = *rState[iFlx];

    // flux for the right and left state (the UPDATE variables have to be passed here, this is implicitly assumed!!!
    m_sumFlux  = updateVarSet->getFlux()(lUpdState, normal[iFlx]);
    m_sumFlux += updateVarSet->getFlux()(rUpdState, normal[iFlx]);

    // compute the Riemann flux
    m_multiRFlux[iFlx] = 0.5*(m_sumFlux - m_absJacob*(rUpdState - lUpdState));

    // increase stateItr
    solStateItr = lastSolState;
  }

  return m_multiRFlux;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& RoeFlux::computeNumDamping(vector< Framework::State* >& lState,
                                                 vector< Framework::State* >& rState,
                                                 const vector< RealVector >& normal,
                                                 const CFuint nbrFlxPnts)
{
  cf_assert(nbrFlxPnts <= lState.size());
  cf_assert(nbrFlxPnts <= rState.size());
  cf_assert(nbrFlxPnts <= normal.size());

  // store the left and right states in vector
  CFuint flxID = 0;
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_updateStates[flxID] = lState[iFlx]; ++flxID;
    m_updateStates[flxID] = rState[iFlx]; ++flxID;
  }

  // get the update VarSet
  SafePtr< ConvectiveVarSet > updateVarSet = getMethodData().getUpdateVar();

  // transform from update states (which are stored) to solution states (in which the equations are written)
 /// @todo broken after release 2009.3
//   m_solStates = m_updateToSolutionVarTrans->transformFromRef(&m_updateStates);

  // compute the numerical dissipation for each flux point
  vector< Framework::State* >::iterator solStateItr = m_solStates->begin();
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    // last position to be copied
    vector< Framework::State* >::iterator lastSolState = solStateItr + 2;

    // copy the (pointers to) left and right states to auxiliary variable
    copy(solStateItr,lastSolState,m_solStatesInFlxPnt.begin());

    // transform from solution states to linear states
    m_linStates = m_solutionToLinearVarTrans->transform(&m_solStatesInFlxPnt);

    // linearize the states (for instance: compute the Roe averaged values) AND set the physical data of the model! This is done inside this function
    m_linearizer->linearize(*m_linStates);

    // set the eigenvectors and eigenvalues of the linearized jacobian
    // USING THE PHYSICAL DATA THAT WAS SET IN linearize();
    updateVarSet->computeEigenValuesVectors(m_rightEv, m_leftEv, m_eValues, normal[iFlx]);

    // set the abs of the  eigen values (the implementation of this
    // function change if there are entropy or carbuncle fixes)
    SetAbsEigenValues();

    // abs of the jacobian
    m_absJacob = m_rightEv*(m_absEvalues*m_leftEv);

    // dereference left and right update states
    State& lUpdState = *lState[iFlx];
    State& rUpdState = *rState[iFlx];

    // compute the numerical dissipation
    m_multiNumDamping[iFlx] = m_absJacob*(rUpdState - lUpdState);

    // increase stateItr
    solStateItr = lastSolState;
  }

  return m_multiNumDamping;
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
  m_solStatesInFlxPnt.resize(2);
  m_updateStatesInFlxPnt.resize(2);
  m_extraVarsInFlxPnt.resize(2);
  m_updateStates.resize(2*m_maxNbrFlxPnts);
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

}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD
