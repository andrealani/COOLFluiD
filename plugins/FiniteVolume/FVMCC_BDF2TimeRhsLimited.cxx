#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_BDF2TimeRhsLimited.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
 
#ifdef CF_HAVE_SINGLE_EXEC
#include "FiniteVolume/MinModTimeLimiter.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_BDF2TimeRhsLimited,
                      CellCenterFVMData,
                      FiniteVolumeModule>
                      fvmccBDF2TimeRhsLimited("BDF2TimeRhsLimited");

//////////////////////////////////////////////////////////////////////////////

void FVMCC_BDF2TimeRhsLimited::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("TimeLimiter","Name of the time limiter.");
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_BDF2TimeRhsLimited::FVMCC_BDF2TimeRhsLimited(const std::string& name) :
  FVMCC_PseudoSteadyTimeRhs(name),
  socket_pastTimeRhs("pastTimeRhs")
{

  addConfigOptionsTo(this);

  _timeLimiterStr = "MinMod";
  setParameter("TimeLimiter",&_timeLimiterStr);

}

//////////////////////////////////////////////////////////////////////////////

FVMCC_BDF2TimeRhsLimited::~FVMCC_BDF2TimeRhsLimited()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_BDF2TimeRhsLimited::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CellCenterFVMCom::configure(args);

#ifndef CF_HAVE_SINGLE_EXEC
  SafePtr<TimeLimiter::PROVIDER> prov = CFNULL;
  try {  
    prov = FACTORY_GET_PROVIDER(getFactoryRegistry(), TimeLimiter, _timeLimiterStr);
  } 
 catch (Common::NoSuchValueException& e) {
   _timeLimiterStr = "Null";
   CFLog(VERBOSE, "FVMCC_BDF2TimeRhsLimited::configure() => choosing NullTimeLimiter instead ..." << "\n");
    prov = FACTORY_GET_PROVIDER(getFactoryRegistry(), TimeLimiter, _timeLimiterStr);
  }
 cf_assert(prov.isNotNull());
 _timeLimiter = prov->create(_timeLimiterStr);
#else
 SelfRegistPtr<TimeLimiter> ptr(new MinModTimeLimiter("MinMod"), CFNULL);
 _timeLimiter .reset(ptr);
#endif 

 cf_assert(_timeLimiter.isNotNull());
 configureNested ( _timeLimiter.getPtr(), args );
}


//////////////////////////////////////////////////////////////////////////////

void FVMCC_BDF2TimeRhsLimited::computeNumericalTransMatrix(const CFuint iState)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> pastTimeRhs = socket_pastTimeRhs.getDataHandle();
  State *const currState = states[iState];

  // this first transformed state HAS TO BE stored,
  // since the returned pointer will change pointee after
  // the second call to transform()
  _tempState = *_updateToSolutionVecTrans->transform(currState);

  const State *const tPastState =
    _updateToSolutionVecTrans->transform(pastStates[iState]);

  const CFreal resFactor = getMethodData().getResFactor();

  // first the contribution to the rhs is computed
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {

    const CFreal dU = _tempState[iEq]*_diagValue - (*tPastState)[iEq]*_diagValue;
    const CFreal slopeNEW = dU;
    const CFreal slopeOLD = pastTimeRhs(iState, iEq, nbEqs);

    CFreal limiterValue = _timeLimiter->computeTimeLimiterValue(slopeOLD,slopeNEW);

    //   CFout <<"TimeLimiter: " << limiterValue <<"\n";
    rhs(iState, iEq, nbEqs) -= dU + limiterValue * ( dU*resFactor - resFactor*pastTimeRhs(iState, iEq, nbEqs));
    if((SubSystemStatusStack::getActive()->isLastStep())
    && (SubSystemStatusStack::getActive()->isSubIterationLastStep()))
    {
      pastTimeRhs(iState, iEq, nbEqs) = dU;
    }
  }

  if ((!getMethodData().isSysMatrixFrozen()) && getMethodData().doComputeJacobian()) {
    // now the contribution to the jacobian matrix is calculated
    _acc->setRowColIndex(0, currState->getLocalID());
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, (*currState)[iVar]);
      
      const RealVector& tempPertState = *_updateToSolutionVecTrans->transform(currState);
      
      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(_tempState,
					 tempPertState,
					 _fluxDiff);
      
      _fluxDiff *= _diagValue*(1.+ resFactor);
      
      _acc->addValues(0, 0, iVar, &_fluxDiff[0]);
      // restore the unperturbed value
      _numericalJacob->restore((*currState)[iVar]);
    }
    
    // add the values in the jacobian matrix
    _lss->getMatrix()->addValues(*_acc);
    
    // reset to zero the entries in the block accumulator
    _acc->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_BDF2TimeRhsLimited::computeAnalyticalTransMatrix(const CFuint iState)
{

  cf_assert(false);

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> pastTimeRhs = socket_pastTimeRhs.getDataHandle();

  State *const currState = states[iState];
  
  const State& tempState =
    *_updateToSolutionVecTrans->transform(currState);
  
  const State& tPastState =
    *_updateToSolutionVecTrans->transform(pastStates[iState]);
  
  const CFreal resFactor = getMethodData().getResFactor();

  // first the contribution to the rhs is computed
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {

  //     const CFreal limiterValue = sqrt(limiter[iState][iEq]);
  //     const CFreal limiterValue = limiter[iState][iEq] * limiter[iState][iEq];
  //    const CFreal limiterValue = limiter[iState][iEq];
  const CFreal limiterValue = 1.;
  const CFreal dU = tempState[iEq] * _diagValue - tPastState[iEq] * _diagValue;

  rhs(iState, iEq, nbEqs) -= dU + limiterValue * ( dU*resFactor - resFactor*pastTimeRhs(iState, iEq, nbEqs));
  if(SubSystemStatusStack::getActive()->isLastStep()) pastTimeRhs(iState, iEq, nbEqs) = dU;
  }
  
  if (getMethodData().doComputeJacobian()) {
    // now the contribution to the jacobian matrix is calculated
    _acc->setRowColIndex(0, currState->getLocalID());
    
    _updateToSolutionInUpdateMatTrans->setMatrix(*states[iState]);
    const RealMatrix& matrix = *_updateToSolutionInUpdateMatTrans->getMatrix();
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      for (CFuint jVar = 0; jVar < nbEqs; ++jVar) {
	const CFreal value = matrix(iVar,jVar)*_diagValue*(1. + resFactor);
	_acc->addValue(0, 0, iVar, jVar, value);
      }
    }
    
    // add the values in the jacobian matrix
    _lss->getMatrix()->addValues(*_acc);
    
    // reset to zero the entries in the block accumulator
    _acc->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
FVMCC_BDF2TimeRhsLimited::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_PseudoSteadyTimeRhs::needsSockets();

  result.push_back(&socket_pastTimeRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
