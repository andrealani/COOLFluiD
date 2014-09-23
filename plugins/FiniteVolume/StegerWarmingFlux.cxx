#include "StegerWarmingFlux.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/BaseTerm.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<StegerWarmingFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
stdStegerWarmingFluxProvider("StegerWarming");

//////////////////////////////////////////////////////////////////////////////

StegerWarmingFlux::StegerWarmingFlux(const std::string& name) :
  FVMCC_FluxSplitter(name),
  _solutionVarSet(CFNULL),
  _updateVarSet(CFNULL),
  _updateToSolutionVarTrans(CFNULL), 
  _solutionStates(CFNULL),
  _statesLR(2),
  _numJacob(new NumericalJacobian("NumericalJacobian")),
  _jacobPlus(),
  _jacobMin(),
  _jacobPlusTrasf(),
  _jacobMinTrasf(),
  _jacobDummy(),
  _eValues(),
  _rightEvalues(),
  _leftEvalues(),
  _tState(),
  _tFlux(),
  _tempUnitNormal()
{
}

//////////////////////////////////////////////////////////////////////////////

StegerWarmingFlux::~StegerWarmingFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void StegerWarmingFlux::compute(RealVector& result)
{
  CellCenterFVMData& data = this->getMethodData(); 
  SafePtr<FVMCC_PolyRec> polyRec = data.getPolyReconstructor();
  
  _statesLR[0] = &polyRec->getCurrLeftState();
  _statesLR[1] = &polyRec->getCurrRightState();
   
  if (!getMethodData().reconstructSolVars()) {
    _solutionStates = getMethodData().getUpdateToSolutionVecTrans()->transform(&_statesLR);
  }
  else {
    _solutionStates = &_statesLR;
  }
  
  // set the jacobian plus and minus
  RealVector& pData = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->getPhysicalData();
  
  _updateVarSet->setExtraData(true);
  

  const RealVector& unitNormal = data.getUnitNormal();
  
  _updateVarSet->computePhysicalData(data.getPolyReconstructor()->getCurrLeftState(), pData);
  _solutionVarSet->splitJacobian(_jacobPlus,
				 _jacobDummy,
				 _eValues,
				 unitNormal);
  
  _updateVarSet->computePhysicalData(data.getPolyReconstructor()->getCurrRightState(), pData);
  _solutionVarSet->splitJacobian(_jacobDummy,
				 _jacobMin,
				 _eValues,
				 unitNormal);
  _updateVarSet->setExtraData(false);
  
  const State& stateL = *(*_solutionStates)[0];
  const State& stateR = *(*_solutionStates)[1];
  // jacobians must be evaluated in different states
  result = _jacobPlus*stateL + _jacobMin*stateR;
  
  // compute update coefficient
  if (!getMethodData().isPerturb()) {
    vector<RealVector>& pdata = polyRec->getExtrapolatedPhysicaData();
    GeometricEntity& face = *data.getCurrentFace();
    const CFreal faceArea = socket_faceAreas.getDataHandle()[face.getID()]/
      polyRec->nbQPoints();
    DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
    
    // left contribution to update coefficient
    CFreal maxEV = _updateVarSet->getMaxEigenValue(pdata[0], unitNormal);
    
    const CFuint leftID = face.getState(0)->getLocalID();
    updateCoeff[leftID] += max(maxEV, (CFreal)0.)*faceArea;
    
    if (!face.getState(1)->isGhost()) {
      // right contribution to update coefficient
      
      _tempUnitNormal = -1.0*unitNormal;
      maxEV = _updateVarSet->getMaxEigenValue(pdata[1],_tempUnitNormal);
      const CFuint rightID = face.getState(1)->getLocalID();
      updateCoeff[rightID] += max(maxEV, (CFreal)0.)*faceArea;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
      
void StegerWarmingFlux::computeLeftJacobian()
{
  CellCenterFVMData& methodData = this->getMethodData();
  State* leftState = methodData.getCurrentFace()->getState(LEFT);
  SafePtr<VarSetMatrixTransformer> vs = methodData.getUpdateToSolutionInUpdateMatTrans();
  vs->setMatrix(*leftState);
  const RealMatrix& dUdP = *vs->getMatrix();
  _lFluxJacobian = _jacobPlus*dUdP;
}
      
//////////////////////////////////////////////////////////////////////////////

void StegerWarmingFlux::computeRightJacobian()
{
  CellCenterFVMData& methodData = this->getMethodData();
  State* rightState = methodData.getCurrentFace()->getState(RIGHT);
  SafePtr<VarSetMatrixTransformer> vs = methodData.getUpdateToSolutionInUpdateMatTrans();
  vs->setMatrix(*rightState);
  const RealMatrix& dUdP = *vs->getMatrix();
  _rFluxJacobian = _jacobMin*dUdP;
}

//////////////////////////////////////////////////////////////////////////////

void StegerWarmingFlux::setup()
{
  FVMCC_FluxSplitter::setup();
  
  _solutionVarSet = getMethodData().getSolutionVar();
  _updateVarSet = getMethodData().getUpdateVar();
  _updateToSolutionVarTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  _jacobPlus.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
		    Framework::PhysicalModelStack::getActive()->getNbEq());
  _jacobMin.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
		   Framework::PhysicalModelStack::getActive()->getNbEq());
  _jacobPlusTrasf.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
			 Framework::PhysicalModelStack::getActive()->getNbEq());
  _jacobMinTrasf.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
			Framework::PhysicalModelStack::getActive()->getNbEq());
  _jacobDummy.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
		     Framework::PhysicalModelStack::getActive()->getNbEq());
  _eValues.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _rightEvalues.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _leftEvalues.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _tState.resize(Framework::PhysicalModelStack::getActive()->getNbEq()); 
  _tFlux.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _tempUnitNormal.resize(Framework::PhysicalModelStack::getActive()->getDim());
  
  RealVector refValues = 
    PhysicalModelStack::getActive()->getImplementor()->getRefStateValues();
  _numJacob->setRefValues(refValues);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
