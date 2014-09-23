#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacobFastCoupling.hh"
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

MethodCommandProvider<FVMCC_ComputeRhsJacobFastCoupling,
		      CellCenterFVMData,
		      FiniteVolumeModule>
fvmcc_computeRhsJacobFastCoupling("NumJacobFastCoupling");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobFastCoupling::FVMCC_ComputeRhsJacobFastCoupling
(const std::string& name) :
  FVMCC_ComputeRhsJacobCoupling(name)
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobFastCoupling::~FVMCC_ComputeRhsJacobFastCoupling()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobFastCoupling::execute()
{
  CFLogDebugMin( "FVMCC_ComputeRhsJacobFastCoupling::execute()" << "\n");

  //reset rhs to 0
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs = 0.0;
  
  // reset to zero all non zero entries in the jacobian
  if (getMethodData().doComputeJacobian()) {
    const CFuint nbLSS = _lss.size();
    for (CFuint i = 0; i < nbLSS; ++i) {
      _lss[i]->getMatrix()->resetToZeroEntries();
    }
  }
  
  _polyRec->computeGradients();
  _polyRec->computeLimiters();

  if (_fluxData->hasDiffusiveTerm) {
    _nodalExtrapolator->extrapolateInAllNodes();
  }

  // set the list of faces
  vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trs.size();

  _faceIdx = 0;

   Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
    geoBuilder = getMethodData().getFaceTrsGeoBuilder();

  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  vector<bool> zeroGrad(PhysicalModelStack::getActive()->getNbEq(), false);

  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];

    CFLog(VERBOSE, "TRS name = " << currTrs->getName() << "\n");

    // the faces on the boundary of the partition don't have to
    // be processed (their fluxes could give NaN)
    if (currTrs->getName() != "PartitionFaces" &&
        currTrs->getName() != "InnerCells") {
      if (currTrs->hasTag("writable")) {
        _currBC = _bcMap.find(iTRS);

        geoData.isBFace = true;
	// set the flags specifying the variables for which the boundary condition
	// imposes constant extrapolation (zero gradient)
	_polyRec->setZeroGradient(_currBC->getZeroGradientsFlags());
      }
      else {
        geoData.isBFace = false;
	_polyRec->setZeroGradient(&zeroGrad);
      }

      // set the current TRS in the geoData
      geoData.trs = currTrs;

      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace, ++_faceIdx) {
        CFLogDebugMed( "iFace = " << iFace << "\n");

	// reset the equation subsystem descriptor
	PhysicalModelStack::getActive()->resetEquationSubSysDescriptor();

        // build the GeometricEntity
        geoData.idx = iFace;
        _currFace = geoBuilder->buildGE();

	// set the data for the FaceIntegrator
	setFaceIntegratorData();

	// extrapolate (and LIMIT, if the reconstruction is linear or more)
	// the solution in the quadrature points
	_polyRec->extrapolate(_currFace);

	computeAndBackUpStatesData();

	const bool isBFace = _currFace->getState(1)->isGhost();
	if (!isBFace) {
	  _fluxSplitter->computeFlux(_flux);
	}
	else {
	  _currBC->computeFlux(_flux);
	}

	if (_fluxData->hasDiffusiveTerm) {
	  // reset to false the flag telling to freeze the diffusive coefficients
	  _diffVar->setFreezeCoeff(false);
	  _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
	  _diffusiveFlux->computeFlux(_dFlux);
	  _flux -= _dFlux;
	}

	//compute the contribution to the RHS
	updateRHS();

	// set on the perturbation flag
	_fluxData->isPerturb = true;

	if (_freezeDiffCoeff) {
	  _diffVar->setFreezeCoeff(true);
	}
	
	if (getMethodData().doComputeJacobian()) {
	  if (!isBFace) {
	    computeJacobianTerm();
	  }
	  else {
	    FVMCC_ComputeRhsJacobCoupling::computeBoundaryJacobianTerm();
	  }
	}

	// set off the perturbation flag
	_fluxData->isPerturb = false;
	geoBuilder->releaseGE();
      }
    }
  }

  // reset the flag for freezing the transport properties
  _diffVar->setFreezeCoeff(false);

  if (getMethodData().isAxisymmetric() || getMethodData().hasSourceTerm()) {
    computeSourceTermContribution();
  }

  // reset the equation subsystem descriptor
  PhysicalModelStack::getActive()->resetEquationSubSysDescriptor();

  // compute dynamically the CFL
  getMethodData().getCFL()->setUpdateOn();
  getMethodData().getCFL()->update();
  getMethodData().getCFL()->setUpdateOff();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobFastCoupling::computeBothJacobTerms(CFuint iLSS)
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);
  BlockAccumulator& acc = *_acc[iLSS];

  // contribution to the jacobian of the first LSS
  // set row and column index for the first block accumulator (LSS 0)
  acc.setRowColIndex(0, state0.getLocalID());
  acc.setRowColIndex(1, state1.getLocalID());
  const CFreal invR0 = _rMid/(state0.getCoordinates())[YY];
  const CFreal invR1 = _rMid/(state1.getCoordinates())[YY];

  SafePtr<vector<CFuint> > equations = _equations[iLSS];
  for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
    // variable to perturb
    const CFuint iVar = (*equations)[iEq];
    CFLogDebugMin("Perturbing iVar = " << iVar << "\n");

    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, state0[iVar]);

    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRec->extrapolate(_currFace, iVar, LEFT);


    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computePerturbedLeftStatesData(iVar);


    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;

    // linearization will be done in the flux splitter if needed
    _fluxSplitter->computeFlux(_pertFlux);

    if (_fluxData->hasDiffusiveTerm) {
      _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);

      // modify only the first _nbEqs components
      _pertFlux.slice(_start) -= _dFlux.slice(_start);
    }

    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(_flux.slice(_start),
				       _pertFlux.slice(_start),
				       _fluxDiff.slice(_start));

    _fluxDiff.slice(_start) *= getResFactor();

    if (!getMethodData().isAxisymmetric()) {
      acc.addValues(0, 0, iEq, &_fluxDiff[_start]);

      // flux is opposite in sign for the other state
      _fluxDiff.slice(_start) *= -1.;

      acc.addValues(1, 0, iEq, &_fluxDiff[_start]);
    }
    else {
      _axiFlux.slice(_start) = _fluxDiff.slice(_start)*invR0;
      acc.addValues(0, 0, iEq, &_axiFlux[_start]);

      _axiFlux.slice(_start) = _fluxDiff.slice(_start)*(-invR1);
      acc.addValues(1, 0, iEq, &_axiFlux[_start]);
    }

    // restore the unperturbed value
    _numericalJacob->restore(state0[iVar]);
    _polyRec->restoreValues(iVar, LEFT);
  }

  restoreLeftStates();

  // second node contribution
  for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
    // variable to perturb
    const CFuint iVar = (*equations)[iEq];
    CFLogDebugMin("Perturbing iVar = " << iVar << "\n");

    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, state1[iVar]);

    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRec->extrapolate(_currFace, iVar, RIGHT);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computePerturbedRightStatesData(iVar);

    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;

    // linearization will be done in the flux splitter if needed
    _fluxSplitter->computeFlux(_pertFlux);

    if (_fluxData->hasDiffusiveTerm) {
      _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);
      _pertFlux.slice(_start) -= _dFlux.slice(_start);
    }

    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(_flux.slice(_start),
				       _pertFlux.slice(_start),
				       _fluxDiff.slice(_start));

    // multiply for the residual factor
    _fluxDiff.slice(_start) *= getResFactor();

    if (!getMethodData().isAxisymmetric()) {
      acc.addValues(0, 1, iEq, &_fluxDiff[_start]);

      // flux is opposite in sign for the other state
      _fluxDiff.slice(_start) *= -1.;

      acc.addValues(1, 1, iEq, &_fluxDiff[_start]);
    }
    else {
      _axiFlux.slice(_start) = _fluxDiff.slice(_start)*invR0;
      acc.addValues(0, 1, iEq, &_axiFlux[_start]);

      _axiFlux.slice(_start) = _fluxDiff.slice(_start)*(-invR1);
      acc.addValues(1, 1, iEq, &_axiFlux[_start]);
    }

    // restore the unperturbed value
    _numericalJacob->restore(state1[iVar]);
    _polyRec->restoreValues(iVar, RIGHT);
  }

  // there is no need to call restoreRightStates();
  // in order to restore the right states because nobody
  // will use them again within this iteration on face
  restoreRightStates();

  // add the values in the jacobian matrix
  _lss[iLSS]->getMatrix()->addValues(acc);

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobFastCoupling::computeJacobTerm(CFuint iLSS,
							 CFuint idx)
{
  cf_assert(_currFace->getState(idx)->isParUpdatable());

  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  // unused // const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal invR   = _rMid/(_currFace->getState(idx)->getCoordinates())[YY];
  const CFreal factor = pow(-1.,static_cast<CFreal>(idx))*getResFactor();

  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);

  BlockAccumulator& acc = *_acc[iLSS];
  acc.setRowColIndex(0, state0.getLocalID());
  acc.setRowColIndex(1, state1.getLocalID());

  // first node contribution
  SafePtr<vector<CFuint> > equations = _equations[iLSS];
  for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
    // variable to perturb
    const CFuint iVar = (*equations)[iEq];
    CFLogDebugMin( "Perturbing iVar = " << iVar << "\n");

    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, state0[iVar]);

    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRec->extrapolate(_currFace, iVar, LEFT);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computePerturbedLeftStatesData(iVar);

    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;

    // linearization will be done in the flux splitter if needed
    _fluxSplitter->computeFlux(_pertFlux);

    if (_fluxData->hasDiffusiveTerm) {
      _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);

      // modify only the first _nbEqs components
      _pertFlux.slice(_start) -= _dFlux.slice(_start);
    }

    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(_flux.slice(_start),
				       _pertFlux.slice(_start),
				       _fluxDiff.slice(_start));

    //multiply for the residual factor
    _fluxDiff.slice(_start) *= factor;

    if (!getMethodData().isAxisymmetric()) {
      acc.addValues(idx, 0, iEq, &_fluxDiff[_start]);
    }
    else {
      _axiFlux.slice(_start) = _fluxDiff.slice(_start)*invR;
      acc.addValues(idx, 0, iEq, &_axiFlux[_start]);
    }

    // restore the unperturbed value
    _numericalJacob->restore(state0[iVar]);
    _polyRec->restoreValues(iVar, LEFT);
  }

  restoreLeftStates();

  // second node contribution
  for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
    // variable to perturb
    const CFuint iVar = (*equations)[iEq];
    CFLogDebugMin( "Perturbing iVar = " << iVar << "\n");

    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, state1[iVar]);
    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points

    _polyRec->extrapolate(_currFace, iVar, RIGHT);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computePerturbedRightStatesData(iVar);

    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;

    // linearization will be done in the flux splitter if needed
    _fluxSplitter->computeFlux(_pertFlux);

    if (_fluxData->hasDiffusiveTerm) {
      _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);

      // modify only the first _nbEqs components
      _pertFlux.slice(_start) -= _dFlux.slice(_start);
    }

    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(_flux.slice(_start),
				       _pertFlux.slice(_start),
				       _fluxDiff.slice(_start));

    // multiply for the residual factor
    _fluxDiff.slice(_start) *= factor;

    if (!getMethodData().isAxisymmetric()) {
      acc.addValues(idx, 1, iEq, &_fluxDiff[_start]);
    }
    else {
      _axiFlux.slice(_start) = _fluxDiff.slice(_start)*invR;
      acc.addValues(idx, 1, iEq, &_axiFlux[_start]);
    }

    // restore the unperturbed value
    _numericalJacob->restore(state1[iVar]);
    _polyRec->restoreValues(iVar, RIGHT);
  }

  // there is no need to call restoreRightStates();
  // in order to restore the right states because nobody
  // will use them again within this iteration on face
  restoreRightStates();

  // add the values in the jacobian matrix
  _lss[iLSS]->getMatrix()->addValues(acc);

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobFastCoupling::computeBoundaryJacobianTerm(CFuint iLSS)
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());

  cf_assert(iLSS < _bAcc.size());

  BlockAccumulator& bAcc = *_bAcc[iLSS];
  SafePtr<LinearSystemSolver> lss = _lss[iLSS];
  const CFreal invR0 = _rMid/(currState.getCoordinates())[YY];

  // we work with a slice [_start, _nbEqs)
  SafePtr<vector<CFuint> > equations = _equations[iLSS];

  if (currState.isParUpdatable()) {
    // copy the original value of the ghost state
    _origState.slice(_start) = ghostState.slice(_start);

    bAcc.setRowColIndex(0, currState.getLocalID());

    for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
      // variable to perturb
      const CFuint iVar = (*equations)[iEq];
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

      if (_fluxData->hasDiffusiveTerm) {
	_nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
        _diffusiveFlux->computeFlux(_dFlux);
        _pertFlux.slice(_start) -= _dFlux.slice(_start);
      }

      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(_flux.slice(_start),
					 _pertFlux.slice(_start),
					 _fluxDiff.slice(_start));

      // multiply for the residual factor
      _fluxDiff.slice(_start) *= getResFactor();

      if (!getMethodData().isAxisymmetric()) {
        bAcc.addValues(0, 0, iEq, &_fluxDiff[_start]);
      }
      else {
        _axiFlux.slice(_start) = _fluxDiff.slice(_start)*invR0;
        bAcc.addValues(0, 0, iEq, &_axiFlux[_start]);
      }

      // restore the unperturbed value
      _numericalJacob->restore(currState[iVar]);

      // restore the original ghost state
      ghostState.slice(_start) = _origState.slice(_start);
    }

    // add the values in the jacobian matrix
    lss->getMatrix()->addValues(bAcc);

    // reset to zero the entries in the block accumulator
    bAcc.reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
