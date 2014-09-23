#include <iterator>

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/CFSide.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/LaxFriedrichsFlux.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    LaxFriedrichsFlux,SpectralFVMethodData,RiemannFlux,SpectralFVModule >
  LaxFriedrichsFluxProvider("LaxFriedrichsFlux");

//////////////////////////////////////////////////////////////////////////////

LaxFriedrichsFlux::LaxFriedrichsFlux(const std::string& name) :
  RiemannFlux(name),
  m_updateToSolutionVarTrans(),
  m_sumFlux(),
  m_updateStates(),
  m_updateStatesInFlxPnt(),
  m_extraVarsInFlxPnt(),
  m_solStates(),
  m_solStatesInFlxPnt()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

LaxFriedrichsFlux::~LaxFriedrichsFlux()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& LaxFriedrichsFlux::computeFlux(State& lState, RealVector& lExtraVars,
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

  // compute the states data associated with lState and rState
  // (must be called before computing the flux and computing the maximum absolute eigenvalue)
 /// @todo broken after release 2009.3
//   updateVarSet->computeStatesData(m_updateStatesInFlxPnt,m_extraVarsInFlxPnt,2);

  // flux for right and left state (the UPDATE variables have to be passed here, this is implicitly assumed!!!
  m_sumFlux  = updateVarSet->getFlux()(lState, normal);
  m_sumFlux += updateVarSet->getFlux()(rState, normal);

  // compute left and right maximum absolute eigenvalues
  const CFreal lMaxAbsEVal = updateVarSet->getMaxAbsEigenValue(lState,normal);
  const CFreal rMaxAbsEVal = updateVarSet->getMaxAbsEigenValue(rState,normal);

  // compute numerical damping coefficient
  const CFreal alphaDmp = 0.5*(lMaxAbsEVal+rMaxAbsEVal);

  // get solution states
  State& lSolState = *(*m_solStates)[LEFT ];
  State& rSolState = *(*m_solStates)[RIGHT];

  // compute the Riemann flux
  m_rFlux = 0.5*(m_sumFlux - alphaDmp*(rSolState - lSolState));

  return m_rFlux;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& LaxFriedrichsFlux::computeNumDamping(State& lState, RealVector& lExtraVars,
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

  // compute the states data associated with lState and rState
  // (must be called before computing the maximum absolute eigenvalue)
 /// @todo broken after release 2009.3
//   updateVarSet->computeStatesData(m_updateStatesInFlxPnt,m_extraVarsInFlxPnt,2);

  // compute left and right maximum absolute eigenvalues
  const CFreal lMaxAbsEVal = updateVarSet->getMaxAbsEigenValue(lState,normal);
  const CFreal rMaxAbsEVal = updateVarSet->getMaxAbsEigenValue(rState,normal);

  // compute numerical damping coefficient
  const CFreal alphaDmp = 0.5*(lMaxAbsEVal+rMaxAbsEVal);

  // get solution states
  State& lSolState = *(*m_solStates)[LEFT ];
  State& rSolState = *(*m_solStates)[RIGHT];

  // compute the numerical damping
  m_numDamping = alphaDmp*(rSolState - lSolState);

  return m_numDamping;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& LaxFriedrichsFlux::computeFlux(vector< Framework::State* >& lState,
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
    // dereference left and right update states
    State& lUpdState = *lState[iFlx];
    State& rUpdState = *rState[iFlx];

    // flux for right and left state (the UPDATE variables have to be passed here, this is implicitly assumed!!!
    m_sumFlux  = updateVarSet->getFlux()(lUpdState,normal[iFlx]);
    m_sumFlux += updateVarSet->getFlux()(rUpdState,normal[iFlx]);

    // compute left and right maximum absolute eigenvalues
    const CFreal lMaxAbsEVal = updateVarSet->getMaxAbsEigenValue(lUpdState,normal[iFlx]);
    const CFreal rMaxAbsEVal = updateVarSet->getMaxAbsEigenValue(rUpdState,normal[iFlx]);

    // compute numerical damping coefficient
    const CFreal alphaDmp = 0.5*(lMaxAbsEVal+rMaxAbsEVal);

    // dereference left and right solution states
    State& lSolState = *(*solStateItr); ++solStateItr;
    State& rSolState = *(*solStateItr); ++solStateItr;

    // compute the Riemann flux
    m_multiRFlux[iFlx] = 0.5*(m_sumFlux - alphaDmp*(rSolState - lSolState));
  }

  return m_multiRFlux;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& LaxFriedrichsFlux::computeNumDamping(vector< Framework::State* >& lState,
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
    // dereference left and right update states
    State& lUpdState = *lState[iFlx];
    State& rUpdState = *rState[iFlx];

    // compute left and right maximum absolute eigenvalues
    const CFreal lMaxAbsEVal = updateVarSet->getMaxAbsEigenValue(lUpdState,normal[iFlx]);
    const CFreal rMaxAbsEVal = updateVarSet->getMaxAbsEigenValue(rUpdState,normal[iFlx]);

    // compute numerical damping coefficient
    const CFreal alphaDmp = 0.5*(lMaxAbsEVal+rMaxAbsEVal);

    // dereference left and right solution states
    State& lSolState = *(*solStateItr); ++solStateItr;
    State& rSolState = *(*solStateItr); ++solStateItr;

    // compute the numerical damping
    m_multiNumDamping[iFlx] = alphaDmp*(rSolState - lSolState);
  }

  return m_multiNumDamping;
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedrichsFlux::setup()
{
  CFAUTOTRACE;

  RiemannFlux::setup();

  // resize variables
  m_solStatesInFlxPnt.resize(2);
  m_updateStatesInFlxPnt.resize(2);
  m_extraVarsInFlxPnt.resize(2);
  m_updateStates.resize(2*m_maxNbrFlxPnts);
  m_sumFlux.resize(m_nbrEqs);

  // get the name of the physical model
  std::string physicsName = PhysicalModelStack::getActive()->getImplementor()->getConvectiveName();

  // get the varset names
  std::string updateVarName = getMethodData().getUpdateVarStr();
  std::string solutionVarName = getMethodData().getSolutionVarStr();

  // create the update to solution variable transformer
  std::string updateToSolutionVarName
    = VarSetTransformer::getProviderName(physicsName,updateVarName,solutionVarName);
  CFLog(INFO, "LaxFriedrichsFlux::setup() => updateToSolutionVarName = " << updateToSolutionVarName << "\n");
  m_updateToSolutionVarTrans =
    Environment::Factory<VarSetTransformer>::getInstance().getProvider
    (updateToSolutionVarName)->create(PhysicalModelStack::getActive()->getImplementor());
  cf_assert(m_updateToSolutionVarTrans.isNotNull());

  //
  m_updateToSolutionVarTrans->setup(2*m_maxNbrFlxPnts);

}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD
