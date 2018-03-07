#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacobAnalytic.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
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

MethodCommandProvider<FVMCC_ComputeRhsJacobAnalytic,
		      CellCenterFVMData,
		      FiniteVolumeModule>
fvmcc_computeRhsJacobAnalytic("NumJacobAnalytic");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobAnalytic::FVMCC_ComputeRhsJacobAnalytic
(const std::string& name) :
  FVMCC_ComputeRhsJacob(name),
  _tmpJacobL(),
  _tmpJacobR()
{
}
    
//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobAnalytic::~FVMCC_ComputeRhsJacobAnalytic()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobAnalytic::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobAnalytic::setup()
{
  FVMCC_ComputeRhsJacob::setup();
  getMethodData().setUseAnalyticalConvJacob(true);
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _tmpJacobL.resize(nbEqs, nbEqs);
  _tmpJacobR.resize(nbEqs, nbEqs);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobAnalytic::computeBothJacobTerms()
{  
  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);
  (!getMethodData().isAxisymmetric()) ? computeNoAxiUpFactors() : computeAxiUpFactors();
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  _acc->setRowColIndex(0, state0.getLocalID());
  _acc->setRowColIndex(1, state1.getLocalID());

  // analytical jacobian of convective fluxes (left and right state)  
  SafePtr<RealMatrix> convJacobL = _fluxSplitter->getLeftFluxJacob();
  SafePtr<RealMatrix> convJacobR = _fluxSplitter->getRightFluxJacob();
  
  if (_hasDiffusiveTerm) {
    // analytical jacobian of diffusive fluxes (left and right state)  
    SafePtr<RealMatrix> diffLeftJacob  = _diffusiveFlux->getLeftFluxJacob();
    SafePtr<RealMatrix> diffRightJacob = _diffusiveFlux->getRightFluxJacob();
    
    // note the opposite sign of the diffusive jacob compared to convective jacob
    *convJacobL -= *diffLeftJacob;
    *convJacobR -= *diffRightJacob;
  }
  
  (!getMethodData().isAxisymmetric()) ? 
    noaxiBothJacobTerms(*convJacobL, *convJacobR) : 
    axiBothJacobTerms(*convJacobL, *convJacobR); 
  
   // _acc->print(); EXIT_AT(1);
  
  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);
  
  // reset to zero the entries in the block accumulator
  _acc->reset();
  
  _sourceJacobOnCell[0] = false;
  _sourceJacobOnCell[1] = false;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobAnalytic::computeJacobTerm(const CFuint idx)
{
  cf_assert (_currFace->getState(idx)->isParUpdatable());
  
  const CFreal factor = pow(-1.,static_cast<CFreal>(idx))*getResFactor();
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  
  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);
  
  _acc->setRowColIndex(0, state0.getLocalID());
  _acc->setRowColIndex(1, state1.getLocalID());
  
  // analytical jacobian of convective fluxes (left and right state)  
  SafePtr<RealMatrix> convJacobL = _fluxSplitter->getLeftFluxJacob();
  SafePtr<RealMatrix> convJacobR = _fluxSplitter->getRightFluxJacob();
  
  if (_hasDiffusiveTerm) {
    // analytical jacobian of diffusive fluxes (left and right state)  
    SafePtr<RealMatrix> diffLeftJacob = _diffusiveFlux->getLeftFluxJacob();
    SafePtr<RealMatrix> diffRightJacob = _diffusiveFlux->getRightFluxJacob();
    
    // note the opposite sign of the diffusive jacob compared to convective jacob
    *convJacobL -= *diffLeftJacob;
    *convJacobR -= *diffRightJacob;
  }
  
  _upStFactor[idx] = (!getMethodData().isAxisymmetric()) ? -getResFactor() :-getResFactor()*_invr[idx];
  
  (!getMethodData().isAxisymmetric()) ? 
    noaxiJacobTerm(idx, *convJacobL, *convJacobR, factor) : 
    axiJacobTerm(idx, *convJacobL, *convJacobR, factor);
  
  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);
  
  // reset to zero the entries in the block accumulator
  _acc->reset(); 
  
  _sourceJacobOnCell[idx] = false;
}

//////////////////////////////////////////////////////////////////////////////
 
void FVMCC_ComputeRhsJacobAnalytic::noaxiJacobTerm(CFuint idx, 
						   RealMatrix& convJacobL, 
						   RealMatrix& convJacobR,
						   CFreal factor)
{
  convJacobL *= factor;
  convJacobR *= factor;
  
  _acc->addValuesM(idx, 0, convJacobL);
  _acc->addValuesM(idx, 1, convJacobR);
  
  if (_sourceJacobOnCell[idx]) {
    addAnalyticSourceTermJacob(idx, _acc.get());
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobAnalytic::axiJacobTerm(CFuint idx, 
						 RealMatrix& convJacobL, 
						 RealMatrix& convJacobR,
						 CFreal factor)
{
  convJacobL *= factor;
  convJacobR *= factor;
  
  const CFreal invR = _rMid*_invr[idx];
  _tmpJacobL = invR*convJacobL;
  _tmpJacobR = invR*convJacobR;
  
  _acc->addValuesM(idx, 0, _tmpJacobL);
  _acc->addValuesM(idx, 1, _tmpJacobR);
  
  if (_sourceJacobOnCell[idx]) {
    addAnalyticSourceTermJacob(idx, _acc.get());
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobAnalytic::noaxiBothJacobTerms(RealMatrix& convJacobL, 
							RealMatrix& convJacobR)
{
  convJacobL *= getResFactor();
  convJacobR *= getResFactor();
  
  _acc->addValuesM(0, 0, convJacobL);
  _acc->addValuesM(0, 1, convJacobR);
  
  // flux is opposite in sign for the other state
  convJacobL *= -1.;
  convJacobR *= -1.;
  
  _acc->addValuesM(1, 0, convJacobL);
  _acc->addValuesM(1, 1, convJacobR);
  
  for (CFuint idx = 0; idx < 2; ++idx) {
    addAnalyticSourceTermJacob(idx, _acc.get());
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobAnalytic::axiBothJacobTerms(RealMatrix& convJacobL, 
						      RealMatrix& convJacobR)
{ 
  convJacobL *= getResFactor();
  convJacobR *= getResFactor();
  
  const CFreal invR0 = _rMid*_invr[0];
  const CFreal invR1 = _rMid*_invr[1];
  
  convJacobL *= invR0;
  convJacobR *= invR0;
  
  _acc->addValuesM(0, 0, convJacobL);
  _acc->addValuesM(0, 1, convJacobR);
  
  const CFreal rRatio = -invR1/invR0;
  convJacobL *= rRatio;
  convJacobR *= rRatio;
  
  _acc->addValuesM(1, 0, convJacobL);
  _acc->addValuesM(1, 1, convJacobR);
  
  for (CFuint idx = 0; idx < 2; ++idx) {
    if (_sourceJacobOnCell[idx]) {
      addAnalyticSourceTermJacob(idx, _acc.get());
    }   
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
