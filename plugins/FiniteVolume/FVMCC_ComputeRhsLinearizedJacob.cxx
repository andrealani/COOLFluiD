#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCC_ComputeRhsLinearizedJacob.hh"
#include "FiniteVolume/FVMCC_BC.hh"

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

MethodCommandProvider<FVMCC_ComputeRhsLinearizedJacob,
		      CellCenterFVMData,
		      FiniteVolumeModule>
FVMCC_ComputeRhsLinearizedJacob("LinearizedAnalyticalJacob");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsLinearizedJacob::FVMCC_ComputeRhsLinearizedJacob(const std::string& name) :
  FVMCC_ComputeRhsJacob_Linearized(name)
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsLinearizedJacob::~FVMCC_ComputeRhsLinearizedJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsLinearizedJacob::computeBothJacobTerms()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);

  _acc->setRowColIndex(0, state0.getLocalID());
  _acc->setRowColIndex(1, state1.getLocalID());

  // extrapolate (and LIMIT, if the reconstruction is linear or more)
  // the solution in the quadrature points
  _polyRec->extrapolate(_currFace);

  // compute the physical data for each left and right reconstructed
  // state and in the left and right cell centers
  computeStatesData2();

  // linearization will be done in the flux splitter if needed
  SafePtr<RealMatrix> leftJacob = _fluxSplitter->getLeftFluxJacob();
  SafePtr<RealMatrix> rightJacob = _fluxSplitter->getRightFluxJacob();

  if (_fluxData->hasDiffusiveTerm) {
    // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
    *leftJacob -= *(_diffusiveFlux->getLeftFluxJacob());
    *rightJacob -= *(_diffusiveFlux->getRightFluxJacob());
  }

  *leftJacob *= getResFactor();
  *rightJacob *= getResFactor();

  if (!getMethodData().isAxisymmetric()) {
    // add jacobian of left diffusive fluxes
    _acc->addValues(0, 0, *leftJacob);
    _acc->addValues(0, 1, *rightJacob);

    // flux is opposite in sign for the other state
    *leftJacob *= -1.;
    *rightJacob *= -1.;

    // add jacobian of left diffusive fluxes
    _acc->addValues(1, 0, *leftJacob);
    _acc->addValues(1, 1, *rightJacob);
  }
  else {
    CFout << "Axi not implemented\n";
    abort();
  }

  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);

  // reset to zero the entries in the block accumulator
  _acc->reset();
 }

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsLinearizedJacob::computeJacobTerm(const CFuint idx)
{
  cf_assert(_currFace->getState(idx)->isParUpdatable());

  const CFreal factor = pow(-1.,static_cast<CFreal>(idx))*getResFactor();

  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);

  _acc->setRowColIndex(0, state0.getLocalID());
  _acc->setRowColIndex(1, state1.getLocalID());

  // extrapolate (and LIMIT, if the reconstruction is linear or more)
  // the solution in the quadrature points
  _polyRec->extrapolate(_currFace);

  // compute the physical data for each left and right reconstructed
  // state and in the left and right cell centers
  computeStatesData2();

  // linearization will be done in the flux splitter if needed
  SafePtr<RealMatrix> leftJacob = _fluxSplitter->getLeftFluxJacob();
  SafePtr<RealMatrix> rightJacob = _fluxSplitter->getRightFluxJacob();

  if (_fluxData->hasDiffusiveTerm) {
    // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
    *leftJacob -= *(_diffusiveFlux->getLeftFluxJacob());
    *rightJacob -= *(_diffusiveFlux->getRightFluxJacob());
  }

  *leftJacob *= factor;
  *rightJacob *= factor;

  if (!getMethodData().isAxisymmetric()) {
    // add jacobian of left diffusive fluxes
    _acc->addValues(idx, 0, *leftJacob);
    _acc->addValues(idx, 1, *rightJacob);
  }
  else {
    CFout << "Axi not implemented\n";
    abort();
  }

  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);

  // reset to zero the entries in the block accumulator
  _acc->reset();
}


//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsLinearizedJacob::computeBoundaryJacobianTerm()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());

  if (currState.isParUpdatable()) {

    // copy the original value of the ghost state
    _origState = ghostState;

    _bAcc->setRowColIndex(0, currState.getLocalID());

    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRec->extrapolate(_currFace);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computeStatesData2();

    // linearization will be done in the flux splitter if needed
    SafePtr<RealMatrix> leftJacob = _fluxSplitter->getLeftFluxJacob();

    if (_fluxData->hasDiffusiveTerm) {
      // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      *leftJacob -= *(_diffusiveFlux->getLeftFluxJacob());
    }

    *leftJacob *= getResFactor();

    if (!getMethodData().isAxisymmetric()) {
      // add jacobian of left diffusive fluxes
      _acc->addValues(0, 0, *leftJacob);
    }
    else {
      CFout << "Axi not implemented\n";
      abort();
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
