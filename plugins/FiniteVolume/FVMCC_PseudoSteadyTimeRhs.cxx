#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCC_PseudoSteadyTimeRhs.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"

#ifdef CF_HAVE_MPI
#include "Common/MPI/MPIStructDef.hh"
#endif

#ifdef CF_HAVE_CUDA
#include "Framework/CudaTimer.hh"
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

MethodCommandProvider<FVMCC_PseudoSteadyTimeRhs, 
		      CellCenterFVMData, 
		      FiniteVolumeModule> 
fvmccPseudoSteadyTimeRhs("PseudoSteadyTimeRhs");

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhs::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >
    ("useGlobalDT", "Flag telling if to use global DT.");
  options.addConfigOption< bool >
    ("useAnalyticalMatrix", "Flag telling if to use analytical matrix.");  
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_PseudoSteadyTimeRhs::FVMCC_PseudoSteadyTimeRhs
(const std::string& name) :
  FVMCC_StdComputeTimeRhs(name),
  _numericalJacob(CFNULL),
  socket_pastStates("pastStates"),
  socket_volumes("volumes"),
  _updateToSolutionVecTrans(CFNULL),
  _updateToSolutionInUpdateMatTrans(CFNULL),
  _fluxDiff(),
  _tempState(),
  _tempPertState(),
  _acc(CFNULL),
  _diagValue(0.0)
{
  addConfigOptionsTo(this);

  _useGlobalDT = false;
  setParameter("useGlobalDT",&_useGlobalDT);

  _useAnalyticalMatrix = false;
  setParameter("useAnalyticalMatrix",&_useAnalyticalMatrix);
}
      
//////////////////////////////////////////////////////////////////////////////

FVMCC_PseudoSteadyTimeRhs::~FVMCC_PseudoSteadyTimeRhs()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhs::setup()
{
  FVMCC_StdComputeTimeRhs::setup();

  _fluxDiff.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempState.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempPertState.resize(PhysicalModelStack::getActive()->getNbEq());

  // numerical jacobian
  _numericalJacob = &getMethodData().getNumericalJacobian();

  _updateToSolutionVecTrans =
    getMethodData().getUpdateToSolutionVecTrans();

  _updateToSolutionInUpdateMatTrans =
    getMethodData().getUpdateToSolutionInUpdateMatTrans();

  _updateToSolutionInUpdateMatTrans->setup(2);
  
  _acc.reset(_lss->createBlockAccumulator(1, 1, PhysicalModelStack::getActive()->getNbEq()));
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhs::execute()
{
  CFAUTOTRACE;

#ifdef CF_HAVE_CUDA
  CudaEnv::CudaTimer& timer = CudaEnv::CudaTimer::getInstance();
  timer.start();
#endif
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbStates = states.size();
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  CFreal minDt = 0.0;
  if (_useGlobalDT) {
    // select the minimum delta T
    minDt = volumes[0]/updateCoeff[0]*cfl;
    for (CFuint iState = 1; iState < nbStates; ++iState) {
      minDt = min(volumes[iState]/updateCoeff[iState]*cfl, minDt);
    }
  }
  
#ifdef CF_HAVE_MPI
  const std::string nsp = getMethodData().getNamespace();
  CFreal totalMinDt = 1e10;
  MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&totalMinDt);
  MPI_Allreduce(&minDt, &totalMinDt, 1, MPI_CFREAL, MPI_MIN, PE::GetPE().GetCommunicator(nsp));
  cf_assert(totalMinDt <= minDt);
  minDt = totalMinDt;
#endif
  
  // if global DT is requested, set the minimum DT
  if (_useGlobalDT) {
    // AL: this is kinda hack, but minDT needs access to CFL, update coefficient and volumes,
    // so it would be more cumbersome to do it elsewhere. Ideally, should be implemented inside a
    // subclass of ComputeDT (future work!!!)
    if (SubSystemStatusStack::getActive()->getSubIter() == 0) {
      SubSystemStatusStack::getActive()->setDT(minDt);
      CFLog(INFO, "FVMCC_PseudoSteadyTimeRhs::execute() => CFL[" << cfl << "], DT[" << minDt <<"]\n");
    }
  }
  
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  
  // add the diagonal entries in the jacobian (updateCoeff/CFL)
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    if (!_useGlobalDT) {
      // dt > 0 for unsteady cases, not for steady!
      _diagValue = (dt > 0.0) ? volumes[iState]/dt : updateCoeff[iState]/cfl;
    }
    else {
      // this corresponds to the case with variable DT, as function of the given CFL 
      cf_assert(minDt > 0.);
      cf_assert(dt > 0.);
      // _diagValue = volumes[iState]/(minDt*cfl); WRONG!!
      _diagValue = volumes[iState]/minDt;      
    }
    
    if (states[iState]->isParUpdatable()) {
      if (!_useAnalyticalMatrix) {
	// compute the transformation matrix numerically
	computeNumericalTransMatrix(iState);
      }
      else {
	computeAnalyticalTransMatrix(iState);
      }
    }
  }

  /* const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  if (nbEqs == 9) {
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  ofstream fout("rhs_after.dat");
  for (CFuint iState = 0; iState < states.size(); ++iState) {
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      fout.precision(14); fout.setf(ios::scientific,ios::floatfield); fout << rhs(iState, iEq, nbEqs) << " ";
    }
    fout << endl;
  }
  }*/
  
#ifdef CF_HAVE_CUDA  
  CFLog(VERBOSE, "FVMCC_PseudoSteadyTimeRhs::execute() took " << timer.elapsed() << " s\n");
#endif
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhs::computeNumericalTransMatrix(const CFuint iState)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();

  State *const currState = states[iState];

  // this first transformed state HAS TO BE stored,
  // since the returned pointer will change pointee after
  // the second call to transform()
  _tempState = static_cast<RealVector&>
    (*_updateToSolutionVecTrans->transform(currState));
  
  const State *const tPastState =
    _updateToSolutionVecTrans->transform(pastStates[iState]);

  //CFLog(INFO, "currState  = " << *currState << "\n");
  //CFLog(INFO, "_tempState = " << _tempState << "\n");
  //CFLog(INFO, "tPastState = " << *tPastState << "\n");
  //CFLog(INFO, "pastStates = " << (*pastStates[iState]) << "\n");
  
  const CFreal resFactor = getMethodData().getResFactor();
  
  // first the contribution to the rhs is computed
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    const CFreal dU = _tempState[iEq] - (*tPastState)[iEq];
    rhs(iState, iEq, nbEqs) -= (!_zeroDiagValue[iEq]) ? dU*_diagValue*resFactor : 0.0;
    /*if (nbEqs == 9) {
      CFLog(INFO, rhs(iState, iEq, nbEqs) << ", dU[" << iEq << "]" << dU << ", resFactor = "
      << resFactor << ", _diagValue = " << _diagValue << "\n");
      }*/
  }
    
  if ((!getMethodData().isSysMatrixFrozen()) && getMethodData().doComputeJacobian()) {
    // now the contribution to the jacobian matrix is calculated
    _acc->setRowColIndex(0, currState->getLocalID());
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, (*currState)[iVar]);
      
      const RealVector& tempPertState = static_cast<RealVector&>
	(*_updateToSolutionVecTrans->transform(currState));
      
      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(_tempState,
					 tempPertState,
					 _fluxDiff);
      
      // _fluxDiff corresponds to a column vector of the dU/dP matrix
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	_fluxDiff[iEq] *= (!_zeroDiagValue[iEq]) ? _diagValue*resFactor : 0.0;
      }
      
      _acc->addValues(0, 0, iVar, &_fluxDiff[0]);
      // restore the unperturbed value
      _numericalJacob->restore((*currState)[iVar]);
    }  
    
    // add the values in the jacobian matrix
    getMethodData().getLSSMatrix(0)->addValues(*_acc);
    
    // reset to zero the entries in the block accumulator
    _acc->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhs::computeAnalyticalTransMatrix
(const CFuint iState)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();

  State *const currState = states[iState];

  const State& tempState =
    *_updateToSolutionVecTrans->transform(currState);

  const State& tPastState =
    *_updateToSolutionVecTrans->transform(pastStates[iState]);

  const CFreal resFactor = getMethodData().getResFactor();

  // first the contribution to the rhs is computed
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    const CFreal dU = tempState[iEq] - tPastState[iEq];
    rhs(iState, iEq, nbEqs) -= (!_zeroDiagValue[iEq]) ? dU*_diagValue*resFactor : 0.0;
  }

  if ((!getMethodData().isSysMatrixFrozen()) && getMethodData().doComputeJacobian()) {
    // now the contribution to the jacobian matrix is calculated
    _acc->setRowColIndex(0, currState->getLocalID());

    _updateToSolutionInUpdateMatTrans->setMatrix(*states[iState]);
    const RealMatrix& matrix = *_updateToSolutionInUpdateMatTrans->getMatrix();
    
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      for (CFuint jVar = 0; jVar < nbEqs; ++jVar) {
	if (!_zeroDiagValue[jVar]) {
	  const CFreal value = matrix(iVar,jVar)*_diagValue*resFactor;
	  _acc->addValue(0, 0, iVar, jVar, value);
	}
      }
    }  
    
    // add the values in the jacobian matrix
    getMethodData().getLSSMatrix(0)->addValues(*_acc);
    
    // reset to zero the entries in the block accumulator
    _acc->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
FVMCC_PseudoSteadyTimeRhs::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_StdComputeTimeRhs::needsSockets();

  result.push_back(&socket_pastStates);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
