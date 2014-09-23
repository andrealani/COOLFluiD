#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacobConv.hh"
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

MethodCommandProvider<FVMCC_ComputeRhsJacobConv,
		      CellCenterFVMData,
		      FiniteVolumeModule>
fvmcc_computeRhsJacobConv("NumJacobConv");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobConv::FVMCC_ComputeRhsJacobConv
(const std::string& name) :
  FVMCC_ComputeRhsJacob(name),
  _tmpJacobL(),
  _tmpJacobR()
{
}
    
//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobConv::~FVMCC_ComputeRhsJacobConv()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobConv::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobConv::setup()
{
  FVMCC_ComputeRhsJacob::setup();
  getMethodData().setUseAnalyticalConvJacob(true);
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _tmpJacobL.resize(nbEqs, nbEqs);
  _tmpJacobR.resize(nbEqs, nbEqs);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobConv::computeBothJacobTerms()
{
  getMethodData().setUseAnalyticalConvJacob(true);
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);

  _acc->setRowColIndex(0, state0.getLocalID());
  _acc->setRowColIndex(1, state1.getLocalID());
  
  // analytical jacobian of convective fluxes (left and right state)  
  SafePtr<RealMatrix> convJacobL = _fluxSplitter->getLeftFluxJacob();
  SafePtr<RealMatrix> convJacobR = _fluxSplitter->getRightFluxJacob();
  // note the opposite sign than convective jacob
  *convJacobL *= getResFactor();
  *convJacobR *= getResFactor();
  
  const CFreal invR0 = _rMid/(state0.getCoordinates())[YY];
  const CFreal invR1 = _rMid/(state1.getCoordinates())[YY];
  
  if (!getMethodData().isAxisymmetric()) {
    _acc->addValues(0, 0, *convJacobL);
    _acc->addValues(0, 1, *convJacobR);
    
    // flux is opposite in sign for the other state
    *convJacobL *= -1.;
    *convJacobR *= -1.;
        
    _acc->addValues(1, 0, *convJacobL);
    _acc->addValues(1, 1, *convJacobR);
  }
  else {
    *convJacobL *= invR0;
    *convJacobR *= invR0;
        
    _acc->addValues(0, 0, *convJacobL);
    _acc->addValues(0, 1, *convJacobR);
    
    const CFreal rRatio = -invR1/invR0;
    *convJacobL *= rRatio;
    *convJacobR *= rRatio;
    
    _acc->addValues(1, 0, *convJacobL);
    _acc->addValues(1, 1, *convJacobR);
  } 
    
  if (_fluxData->hasDiffusiveTerm) {
    // first node contribution
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
      
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, state0[iVar]);
      
      // set the perturbed variable
      _fluxData->iPerturbVar = iVar;
      
      // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_pertFlux);
                  
      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(_dFlux,_pertFlux,_fluxDiff);
            
      //multiply for the residual factor
      _fluxDiff *= (-getResFactor());
            
      if (!getMethodData().isAxisymmetric()) {
	_acc->addValues(0, 0, iVar, &_fluxDiff[0]);
	
	// flux is opposite in sign for the other state
	_fluxDiff *= -1.;
	
	_acc->addValues(1, 0, iVar, &_fluxDiff[0]);
      }
      else {
	_axiFlux = _fluxDiff*invR0;
	_acc->addValues(0, 0, iVar, &_axiFlux[0]);
	
	_axiFlux = _fluxDiff*(-invR1);
	_acc->addValues(1, 0, iVar, &_axiFlux[0]);
      }
      
      // restore the unperturbed value
      _numericalJacob->restore(state0[iVar]);
      
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, state1[iVar]);
      
      // set the perturbed variable
      _fluxData->iPerturbVar = iVar;
      
      // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_pertFlux);
            	    
      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(_dFlux,_pertFlux,_fluxDiff);
      // multiply for the residual factor
      _fluxDiff *= (-getResFactor());
            
      if (!getMethodData().isAxisymmetric()) {
	_acc->addValues(0, 1, iVar, &_fluxDiff[0]);
	
	// flux is opposite in sign for the other state
	_fluxDiff *= -1.;
	
	_acc->addValues(1, 1, iVar, &_fluxDiff[0]);
      }
      else {
	_axiFlux = _fluxDiff*invR0;
	_acc->addValues(0, 1, iVar, &_axiFlux[0]);
	
	_axiFlux = _fluxDiff*(-invR1);
	_acc->addValues(1, 1, iVar, &_axiFlux[0]);
      }
      
      // restore the unperturbed value
      _numericalJacob->restore(state1[iVar]);
    }
  }
  
  // _acc->print();
  
  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);
  
  // reset to zero the entries in the block accumulator
  _acc->reset();

}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobConv::computeJacobTerm(const CFuint idx)
{
  cf_assert (_currFace->getState(idx)->isParUpdatable());

  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal invR   = _rMid/(_currFace->getState(idx)->getCoordinates())[YY];
  const CFreal factor = pow(-1.,static_cast<CFreal>(idx))*getResFactor();

  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);

  _acc->setRowColIndex(0, state0.getLocalID());
  _acc->setRowColIndex(1, state1.getLocalID());
  
  // analytical jacobian of convective fluxes (left and right state)  
  SafePtr<RealMatrix> convJacobL = _fluxSplitter->getLeftFluxJacob();
  SafePtr<RealMatrix> convJacobR = _fluxSplitter->getRightFluxJacob();
  // note the opposite sign than convective jacob
  *convJacobL *= factor;
  *convJacobR *= factor;
  
  if (!getMethodData().isAxisymmetric()) {
    _acc->addValues(idx, 0, *convJacobL);
    _acc->addValues(idx, 1, *convJacobR);
  }
  else {
    _tmpJacobL = invR*(*convJacobL);
    _tmpJacobR = invR*(*convJacobR);
    
    _acc->addValues(idx, 0, _tmpJacobL);
    _acc->addValues(idx, 1, _tmpJacobR);
  } 
  
  if (_fluxData->hasDiffusiveTerm) {
    // first node contribution
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
      
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, state0[iVar]);
      
      // set the perturbed variable
      _fluxData->iPerturbVar = iVar;
      
      // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_pertFlux);
      
      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(_dFlux,_pertFlux,_fluxDiff);
      
      //multiply for the residual factor
      _fluxDiff *= (-factor);
      
      if (!getMethodData().isAxisymmetric()) {
	_acc->addValues(idx, 0, iVar, &_fluxDiff[0]);
      }
      else {
	_axiFlux = _fluxDiff*invR;
	_acc->addValues(idx, 0, iVar, &_axiFlux[0]);
      }
      
      // restore the unperturbed value
      _numericalJacob->restore(state0[iVar]);
      
      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, state1[iVar]);
      
      // set the perturbed variable
      _fluxData->iPerturbVar = iVar;
      
      // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_pertFlux);
      
      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(_dFlux, _pertFlux, _fluxDiff);
      
      // multiply for the residual factor
      _fluxDiff *= (-factor);
      
      if (!getMethodData().isAxisymmetric()) {
	_acc->addValues(idx, 1, iVar, &_fluxDiff[0]);
      }
      else {
	_axiFlux = _fluxDiff*invR;
	_acc->addValues(idx, 1, iVar, &_axiFlux[0]);
      }
      
      // restore the unperturbed value
      _numericalJacob->restore(state1[iVar]);
    }
  }
  
  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);
  
  // reset to zero the entries in the block accumulator
  _acc->reset();
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
