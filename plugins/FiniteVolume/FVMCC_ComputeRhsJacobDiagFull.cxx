#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacobDiagFull.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_ComputeRhsJacobDiagFull,
		      CellCenterFVMData,
		      FiniteVolumeModule>
fvmcc_computeRhsJacobDiagFull("NumJacobDiagFull");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobDiagFull::FVMCC_ComputeRhsJacobDiagFull
(const std::string& name) :
  FVMCC_ComputeRhsJacobDiag(name),
  _tmpJacob(),
  _avState(CFNULL)
{
}
    
//////////////////////////////////////////////////////////////////////////////
      
FVMCC_ComputeRhsJacobDiagFull::~FVMCC_ComputeRhsJacobDiagFull()
{
  deletePtr(_avState);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiagFull::setup()
{
  FVMCC_ComputeRhsJacobDiag::setup();
  
  getMethodData().setUseAnalyticalConvJacob(true);
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _tmpJacob.resize(nbEqs, nbEqs);
  _avState = new State();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiagFull::computeBothJacobTerms()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;
  
  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);
  
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();
  
  const CFint stateID0 = upLocalIDsAll[state0.getLocalID()];
  const CFint stateID1 = upLocalIDsAll[state1.getLocalID()];
  cf_assert(stateID0 >= 0);
  cf_assert(stateID1 >= 0);
  
  _matIter0.resetPtr(&diagMatrices[stateID0*nbEqs2]);
  _matIter1.resetPtr(&diagMatrices[stateID1*nbEqs2]);
  
  getMethodData().setUseAnalyticalConvJacob(true);
  
  // analytical jacobian of convective fluxes (left and right state)  
  SafePtr<RealMatrix> convJacobL = _fluxSplitter->getLeftFluxJacob();
  SafePtr<RealMatrix> convJacobR = _fluxSplitter->getRightFluxJacob();
  // note the opposite sign than convective jacob
  *convJacobL *= getResFactor();
  *convJacobR *= getResFactor();
  
  if (!getMethodData().isAxisymmetric()) {
    _matIter0 += *convJacobL;
    
    // flux is opposite in sign for the other state
    *convJacobL *= -1.;
    *convJacobR *= -1.;
    
    _matIter1 += *convJacobR;
  }
  else {
    const CFreal invR0 = _rMid/(state0.getCoordinates())[YY];
    const CFreal invR1 = _rMid/(state1.getCoordinates())[YY];
    
    *convJacobL *= invR0;
    *convJacobR *= invR0;
    
    _matIter0  += *convJacobL;
    
    const CFreal rRatio = -invR1/invR0;
    *convJacobL *= rRatio;
    *convJacobR *= rRatio;
    
    _matIter1  += *convJacobR;
  } 
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiagFull::computeJacobTerm(const CFuint idx)
{
  cf_assert(_currFace->getState(idx)->isParUpdatable());
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;
  
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();
  
  RealVector& pData = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->getPhysicalData();
  SafePtr<ConvectiveVarSet> solutionVars = getMethodData().getSolutionVar(); 
  SafePtr<VarSetMatrixTransformer> vs = getMethodData().getUpdateToSolutionInUpdateMatTrans();
  const CFreal factor = 0.5*pow(-1.,static_cast<CFreal>(idx))*getResFactor()*
    (*_fluxData->faceAreas)[_faceIdx];
  
  if (idx == 0) {
    // dFdU (U_0)
    State& state0 = *_currFace->getState(0);
    const CFint stateID0 = upLocalIDsAll[state0.getLocalID()];
    cf_assert(stateID0 >= 0);

    _matIter0.resetPtr(&diagMatrices[stateID0*nbEqs2]);
    
    _reconstrVar->computePhysicalData(state0, pData);
    solutionVars->computeProjectedJacobian(_fluxData->unitNormal, _tmpJacob);
    vs->setMatrix(state0);
    
    const CFreal absLambda0 = std::abs(_reconstrVar->getMaxEigenValue(state0, _fluxData->unitNormal));
    for (CFuint i = 0; i < nbEqs; ++i) {
      _tmpJacob(i,i) += absLambda0;
    }
    _matIter0 = _tmpJacob*(*vs->getMatrix());
    _matIter0 *= factor;
  }
  
  if (idx == 1) { 
    // dFdU (U_1)
    State& state1 = *_currFace->getState(1);
    const CFint stateID1 = upLocalIDsAll[state1.getLocalID()];
    cf_assert(stateID1 >= 0);
    
    _matIter1.resetPtr(&diagMatrices[stateID1*nbEqs2]);
    
    _reconstrVar->computePhysicalData(state1, pData);
    solutionVars->computeProjectedJacobian(_fluxData->unitNormal, _tmpJacob);
    vs->setMatrix(state1);
    const CFreal absLambda1 = std::abs(_reconstrVar->getMaxEigenValue(state1, _fluxData->unitNormal));
    for (CFuint i = 0; i < nbEqs; ++i) {
      _tmpJacob(i,i) += absLambda1;
    }
    
    _matIter1 = _tmpJacob*(*vs->getMatrix());
    _matIter1 *= factor;
  }
}


//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiagFull::computeBoundaryJacobianTerm()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqs2 = nbEqs*nbEqs;
  
  State& state0 = *_currFace->getState(0);
  
  DataHandle<CFreal> diagMatrices = socket_diagMatrices.getDataHandle();
  DataHandle<CFint> upLocalIDsAll = socket_upLocalIDsAll.getDataHandle();

  const CFint stateID0 = upLocalIDsAll[state0.getLocalID()];
  cf_assert(stateID0 >= 0);
  _matIter0.resetPtr(&diagMatrices[stateID0*nbEqs2]);
  
  getMethodData().setUseAnalyticalConvJacob(true);
  
  // analytical jacobian of convective fluxes (left and right state)  
  SafePtr<RealMatrix> convJacobL = _fluxSplitter->getLeftFluxJacob();
  // note the opposite sign than convective jacob
  *convJacobL *= getResFactor();
  
  if (!getMethodData().isAxisymmetric()) {
    _matIter0 += *convJacobL;
  }
  else {
    const CFreal invR0 = _rMid/(state0.getCoordinates())[YY];
    *convJacobL *= invR0;
    _matIter0 += *convJacobL;
  } 
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
