#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCC_PseudoSteadyTimeRhsBlockDiag.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulatorBase.hh"
#include "Framework/MeshData.hh"
#include "MathTools/MatrixInverter.hh"

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
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_PseudoSteadyTimeRhsBlockDiag, 
		      CellCenterFVMData, 
		      FiniteVolumeModule> 
fvmccPseudoSteadyTimeRhsBlockDiag("PseudoSteadyTimeRhsBlockDiag");

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsBlockDiag::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >
    ("useGlobalDT", "Flag telling if to use global DT.");
  options.addConfigOption< bool >
    ("useAnalyticalMatrix", "Flag telling if to use analytical matrix.");  
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_PseudoSteadyTimeRhsBlockDiag::FVMCC_PseudoSteadyTimeRhsBlockDiag
(const std::string& name) :
  FVMCC_StdComputeTimeRhs(name),
  _numericalJacob(CFNULL),
  socket_pastStates("pastStates"),
  socket_volumes("volumes"),
  socket_diagMatrix("diagMatrix"),
  _updateToSolutionVecTrans(CFNULL),
  _updateToSolutionInUpdateMatTrans(CFNULL),
  _inverter(CFNULL),
  _fluxDiff(),
  _tempState(),
  _tempPertState(),
  _tempRHS(),
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

FVMCC_PseudoSteadyTimeRhsBlockDiag::~FVMCC_PseudoSteadyTimeRhsBlockDiag()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsBlockDiag::setup()
{
  FVMCC_StdComputeTimeRhs::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _fluxDiff.resize(nbEqs);
  _tempState.resize(nbEqs);
  _tempPertState.resize(nbEqs);
  _tempRHS.resize(nbEqs*socket_rhs.getDataHandle().size());
  
  // numerical jacobian
  _numericalJacob = &getMethodData().getNumericalJacobian();

  _updateToSolutionVecTrans =
    getMethodData().getUpdateToSolutionVecTrans();

  _updateToSolutionInUpdateMatTrans =
    getMethodData().getUpdateToSolutionInUpdateMatTrans();

  _updateToSolutionInUpdateMatTrans->setup(2);

  _inverter.reset(MatrixInverter::create(nbEqs, false));
  
  _acc.reset(new BlockAccumulatorBase(1, 1, nbEqs, CFNULL));
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsBlockDiag::execute()
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
      CFLog(INFO, "FVMCC_PseudoSteadyTimeRhsBlockDiag::execute() => CFL[" << cfl << "], DT[" << minDt <<"]\n");
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

  solveSys();
  
#ifdef CF_HAVE_CUDA  
  CFLog(VERBOSE, "FVMCC_PseudoSteadyTimeRhsBlockDiag::execute() took " << timer.elapsed() << " s\n");
#endif
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsBlockDiag::computeNumericalTransMatrix(const CFuint iState)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle< CFreal> diagMatrix = socket_diagMatrix.getDataHandle();
  
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
    const CFuint nbEqsSq = nbEqs*nbEqs;
    _acc->reset(1, 1, nbEqs, &diagMatrix[iState*nbEqsSq]);
    
    // now the contribution to the jacobian matrix is calculated
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
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsBlockDiag::computeAnalyticalTransMatrix
(const CFuint iState)
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle< CFreal> diagMatrix = socket_diagMatrix.getDataHandle();
  
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
    const CFuint nbEqsSq = nbEqs*nbEqs;
    _acc->reset(1, 1, nbEqs, &diagMatrix[iState*nbEqsSq]);
    
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
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
FVMCC_PseudoSteadyTimeRhsBlockDiag::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_StdComputeTimeRhs::needsSockets();

  result.push_back(&socket_pastStates);
  result.push_back(&socket_volumes);
  result.push_back(&socket_diagMatrix);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsBlockDiag::solveSys()
{
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> diagMatrix = socket_diagMatrix.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;
  const CFuint nbStates = socket_states.getDataHandle().size();
  
  RealMatrix inverse(nbEqs, nbEqs);
  RealVector b(nbEqs);
  RealMatrix mat(nbEqs, nbEqs, static_cast<CFreal*>(CFNULL));
  RealVector dU(nbEqs, static_cast<CFreal*>(CFNULL));
  
  for (CFuint i = 0; i < rhs.size(); ++i) {_tempRHS[i] = rhs[i];}
  
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    if (states[iState]->isParUpdatable()) {
      const CFuint startv = iState*nbEqs;
      const CFuint startm = iState*nbEqs2;
      
      // assign to the rhs the result: inverse*rhs
      dU.wrap(nbEqs, &rhs[startv]);
      mat.wrap(nbEqs, nbEqs, &diagMatrix[startm]);
      b = dU; // back up rhs
      inverse = 0.;
      _inverter->invert(mat, inverse);
      dU = inverse*b;
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
