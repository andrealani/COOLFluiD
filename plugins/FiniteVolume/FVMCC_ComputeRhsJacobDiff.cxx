#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacobDiff.hh"
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

MethodCommandProvider<FVMCC_ComputeRhsJacobDiff,
		      CellCenterFVMData,
		      FiniteVolumeModule>
fvmcc_computeRhsJacobDiff("NumJacobDiff");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobDiff::FVMCC_ComputeRhsJacobDiff(const std::string& name) :
  FVMCC_ComputeRhsJacob(name)
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobDiff::~FVMCC_ComputeRhsJacobDiff()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiff::computeBothJacobTerms()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);
  
  _acc->setRowColIndex(0, state0.getLocalID());
  _acc->setRowColIndex(1, state1.getLocalID());
  
  // analytical jacobian of diffusive fluxes (left and right state)  
  SafePtr<RealMatrix> diffLeftJacob = _diffusiveFlux->getLeftFluxJacob();
  SafePtr<RealMatrix> diffRightJacob = _diffusiveFlux->getRightFluxJacob();
  // note the opposite sign than convective jacob
  *diffLeftJacob *= -getResFactor();
  *diffRightJacob *= -getResFactor();
  
  if (!getMethodData().isAxisymmetric()) {
    // add jacobian of left diffusive fluxes
    _acc->addValues(0, 0, *diffLeftJacob);
    _acc->addValues(0, 1, *diffRightJacob);
    
    // flux is opposite in sign for the other state
    *diffLeftJacob *= -1.;
    *diffRightJacob *= -1.;
    
    // add jacobian of left diffusive fluxes
    _acc->addValues(1, 0, *diffLeftJacob);
    _acc->addValues(1, 1, *diffRightJacob);
  }
  else {
    cout << "axy not implemented"<< endl; abort();
    //   _axiFlux = _fluxDiff*invR0;
    //     _acc->addValues(0, 0, iVar, &_axiFlux[0]);
    
    //     _axiFlux = _fluxDiff*(-invR1);
    //     _acc->addValues(1, 0, iVar, &_axiFlux[0]);
  } 
  
  const CFreal invR0 = _rMid/(state0.getCoordinates())[YY];
  const CFreal invR1 = _rMid/(state1.getCoordinates())[YY];

  // first node contribution
  for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");

    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, state0[iVar]);

    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRec->extrapolate(_currFace, iVar, LEFT);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computeStatesData();

    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;
    
    _fluxSplitter->computeFlux(_pertFlux);
    
    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(_flux,_pertFlux,_fluxDiff);
    _fluxDiff *= getResFactor();
    
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
    _polyRec->restoreValues(iVar, LEFT);
    
    // perturb the given component of the second state vector
    _numericalJacob->perturb(iVar, state1[iVar]);

    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRec->extrapolate(_currFace, iVar, RIGHT);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computeStatesData();

    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;

    _fluxSplitter->computeFlux(_pertFlux);
    
    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(_flux, _pertFlux, _fluxDiff);
    _fluxDiff *= getResFactor();
    
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
    _polyRec->restoreValues(iVar, RIGHT);
  }
    
  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);
  
  // reset to zero the entries in the block accumulator
  _acc->reset();
 }

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobDiff::computeBoundaryJacobianTerm()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());
  const CFreal invR0 = _rMid/(currState.getCoordinates())[YY];
  
  // analytical jacobian of diffusive fluxes (left state)  
  SafePtr<RealMatrix> diffLeftJacob = _diffusiveFlux->getLeftFluxJacob();
  // note the opposite sign than convective jacob
  *diffLeftJacob *= -getResFactor();
  
  if (!getMethodData().isAxisymmetric()) {
    // add jacobian of left diffusive fluxes
    _bAcc->addValues(0,0, *diffLeftJacob);
  }
  else {
    *diffLeftJacob *= invR0;
    _bAcc->addValues(0, 0, *diffLeftJacob);
  }     
    
  if (currState.isParUpdatable()) {

    // copy the original value of the ghost state
    _origState = ghostState;

    _bAcc->setRowColIndex(0, currState.getLocalID());

    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");

      // perturb the given component of the state vector
      _numericalJacob->perturb(iVar, currState[iVar]);

      // compute the ghost state in the perturbed inner state
      _currBC->setGhostState(_currFace);

      // extrapolate (and LIMIT, if the reconstruction is linear or more)
      // the solution in the quadrature points
      _polyRec->extrapolate(_currFace);

      // compute the physical data for each left and right reconstructed
      // state and in the left and right cell centers
      computeStatesData();

      // set the perturbed variable
      _fluxData->iPerturbVar = iVar;

      _currBC->computeFlux(_pertFlux);

      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(_flux, _pertFlux, _fluxDiff);
      _fluxDiff *= getResFactor();
      
      if (!getMethodData().isAxisymmetric()) {
        _bAcc->addValues(0, 0, iVar, &_fluxDiff[0]);
      }
      else {
	_axiFlux = _fluxDiff*invR0;
        _bAcc->addValues(0, 0, iVar, &_axiFlux[0]);
      }     
      
      // restore the unperturbed value
      _numericalJacob->restore(currState[iVar]);

      // restore the original ghost state
      ghostState = _origState;
    }
    
    // add the values in the jacobian matrix
    _lss->getMatrix()->addValues(*_bAcc);
    
    // reset to zero the entries in the block accumulator
    _bAcc->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
