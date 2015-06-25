#include "RoeFlux.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/BaseTerm.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RoeFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
stdRoeFluxProvider("Roe");
      
//////////////////////////////////////////////////////////////////////////////

void RoeFlux::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal,Config::DynamicOption<> >
    ("DiffCoeff", "Diffusion reduction coefficient");
}

//////////////////////////////////////////////////////////////////////////////

RoeFlux::RoeFlux(const std::string& name) :
  FVMCC_FluxSplitter(name),
  _sumFlux(),
  _rightEv(),
  _leftEv(),
  _eValues(),
  _rightEvalues(),
  _leftEvalues(),
  _absEvalues(),
  _tState(),
  _absJacob(),
  _jRight(),
  _jLeft(),
  _jacob(),
  _jacobLeftTransf(),
  _jacobRightTransf(),
  _jacobDummy(),
  _tempUnitNormal(),
  _solutionStates(CFNULL),
  _statesLR(2)
{
  addConfigOptionsTo(this);
  _currentDiffRedCoeff = 1.0;
  setParameter("DiffCoeff", &_currentDiffRedCoeff);
}
      
//////////////////////////////////////////////////////////////////////////////

RoeFlux::~RoeFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void RoeFlux::configure ( Config::ConfigArgs& args )
{
  FVMCC_FluxSplitter::configure(args);
  CFLog(VERBOSE, "RoeFlux::configure() => DiffCoeff = " << getReductionCoeff() << " \n");
}

//////////////////////////////////////////////////////////////////////////////
 
void RoeFlux::setup()
{
  FVMCC_FluxSplitter::setup();
  
  _sumFlux.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _rightEv.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
		  Framework::PhysicalModelStack::getActive()->getNbEq());
  _leftEv.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
		 Framework::PhysicalModelStack::getActive()->getNbEq());
  _eValues.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _rightEvalues.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _leftEvalues.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _absEvalues.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _tState.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _absJacob.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
		   Framework::PhysicalModelStack::getActive()->getNbEq());
  _jRight.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
		   Framework::PhysicalModelStack::getActive()->getNbEq());
  _jLeft.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
		Framework::PhysicalModelStack::getActive()->getNbEq());
  _jacob.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
		Framework::PhysicalModelStack::getActive()->getNbEq());
  _jacobLeftTransf.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
			  Framework::PhysicalModelStack::getActive()->getNbEq());
  _jacobRightTransf.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
		    Framework::PhysicalModelStack::getActive()->getNbEq());
  
  _jacobDummy.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
		     Framework::PhysicalModelStack::getActive()->getNbEq());
  _tempUnitNormal.resize(Framework::PhysicalModelStack::getActive()->getDim());
  
  RealVector refValues = 
    PhysicalModelStack::getActive()->getImplementor()->getRefStateValues();
  getMethodData().getNumericalJacobian().setRefValues(refValues);
}
      
//////////////////////////////////////////////////////////////////////////////
 
void RoeFlux::compute(RealVector& result)
{  
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  SafePtr<ConvectiveVarSet> solutionVarSet = getMethodData().getSolutionVar();
  CellCenterFVMData& data = this->getMethodData(); 
  GeometricEntity& face = *data.getCurrentFace();
  SafePtr<FVMCC_PolyRec> polyRec = data.getPolyReconstructor();
  
  _statesLR[0] = &polyRec->getCurrLeftState();
  _statesLR[1] = &polyRec->getCurrRightState();
  
  cf_assert(*_statesLR[0] == polyRec->getCurrLeftState());
  cf_assert(*_statesLR[1] == polyRec->getCurrRightState());
  
  if (!getMethodData().reconstructSolVars()) {
    _solutionStates = getMethodData().getUpdateToSolutionVecTrans()->transform(&_statesLR);
  }
  else {
    _solutionStates = &_statesLR;
  }
  
  linearize();
  
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  
  // set the eigenvectors and eigenvalues of the linearized jacobian
  solutionVarSet->computeEigenValuesVectors(_rightEv,
					    _leftEv,
					    _eValues,
					    unitNormal);
  
  // set the abs of the  eigen values (the implementation of this
  // function change if there are entropy or carbuncle fixes)
  setAbsEigenValues();
  
  // flux for the right and left state
  vector<RealVector>& pdata = polyRec->getExtrapolatedPhysicaData();
  _sumFlux =  updateVarSet->getFlux()(pdata[0], unitNormal);
  _sumFlux += updateVarSet->getFlux()(pdata[1], unitNormal);
  
  const State& stateL = *(*_solutionStates)[0];
  const State& stateR = *(*_solutionStates)[1];
  result = 0.5*(_sumFlux - getReductionCoeff()*(_rightEv*(_absEvalues*_leftEv))*(stateR - stateL));
  
  // compute update coefficient
  if (!getMethodData().isPerturb()) {    
    DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
    const CFreal faceArea = socket_faceAreas.getDataHandle()[face.getID()]/
      polyRec->nbQPoints();
    
    // left contribution to update coefficient
    CFreal maxEV = updateVarSet->getMaxEigenValue(pdata[0], unitNormal);
    
    const CFuint leftID = face.getState(0)->getLocalID();
    updateCoeff[leftID] += max(maxEV, (CFreal)0.)*faceArea;
    
    if (!face.getState(1)->isGhost()) {
      // right contribution to update coefficient
      
      _tempUnitNormal = -1.0*unitNormal;
      maxEV = updateVarSet->getMaxEigenValue(pdata[1],_tempUnitNormal);
      
      const CFuint rightID = face.getState(1)->getLocalID();
      updateCoeff[rightID] += max(maxEV, (CFreal)0.)*faceArea;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RoeFlux::linearize()
{
  vector<State*> *const linearStates = getMethodData().getSolutionToLinearVecTrans()->
    transform(_solutionStates);
  getMethodData().getJacobianLinearizer()->linearize(*linearStates);
}
      
//////////////////////////////////////////////////////////////////////////////

void RoeFlux::setAbsEigenValues()
{
  _absEvalues = abs(_eValues);
}

//////////////////////////////////////////////////////////////////////////////

void RoeFlux::computeTransformMatrix(State* currState)
{
  // this first transformed state HAS TO BE stored,
  // since the returned pointer will change pointee after
  // the second call to transform()
  _tState = static_cast<RealVector&>
    (*getMethodData().getUpdateToSolutionVecTrans()->transform(currState));
  
  NumericalJacobian& numJacob = getMethodData().getNumericalJacobian();
  // first the contribution to the rhs is computed
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    // perturb the given component of the state vector
    numJacob.perturb(iVar, (*currState)[iVar]);
    
    const RealVector& tPertState = static_cast<RealVector&>
      (*getMethodData().getUpdateToSolutionVecTrans()->transform(currState));
    
    // compute the finite difference derivative of the flux
    numJacob.computeDerivative(_tState, tPertState, _sumFlux);
    _jacobDummy.setColumn(_sumFlux,iVar);
    
    // restore the unperturbed value
    numJacob.restore((*currState)[iVar]);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void RoeFlux::computeLeftJacobian()
{ 
  CellCenterFVMData& data = this->getMethodData(); 
  GeometricEntity& face = *data.getCurrentFace();
  State *const leftState = face.getState(LEFT);
  RealVector& pData = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->getPhysicalData();
  getMethodData().getUpdateVar()->computePhysicalData(*leftState, pData);
  getMethodData().getSolutionVar()->computeProjectedJacobian(data.getUnitNormal(), _jacob); 
  _jLeft = 0.5*(_jacob + _absJacob);
  
  // computeTransformMatrix(leftState);
  // _lFluxJacobian = _jLeft*_jacobDummy;
  
  SafePtr<VarSetMatrixTransformer>  vs = 
    getMethodData().getUpdateToSolutionInUpdateMatTrans();
  vs->setMatrix(*leftState);
  const RealMatrix& dUdP = *vs->getMatrix();  
  _lFluxJacobian = _jLeft*dUdP;
}
      
//////////////////////////////////////////////////////////////////////////////

void RoeFlux::computeRightJacobian()
{ 
  CellCenterFVMData& data = this->getMethodData(); 
  GeometricEntity& face = *data.getCurrentFace();
  State *const rightState = face.getState(RIGHT);
  RealVector& pData = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->getPhysicalData();
  getMethodData().getUpdateVar()->computePhysicalData(*rightState, pData);
  getMethodData().getSolutionVar()->computeProjectedJacobian(data.getUnitNormal(), _jacob);
  _jRight = 0.5*(_jacob - _absJacob);
  
  // computeTransformMatrix(rightState);
  //_rFluxJacobian = _jRight*_jacobDummy;
  
  SafePtr<VarSetMatrixTransformer>  vs = 
    getMethodData().getUpdateToSolutionInUpdateMatTrans();
  vs->setMatrix(*rightState);
  const RealMatrix& dUdP = *vs->getMatrix();
  _rFluxJacobian = _jRight*dUdP;
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
