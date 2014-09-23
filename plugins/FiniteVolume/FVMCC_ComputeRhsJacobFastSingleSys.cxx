#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacobFastSingleSys.hh"
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

MethodCommandProvider<FVMCC_ComputeRhsJacobFastSingleSys,
		      CellCenterFVMData,
		      FiniteVolumeModule>
FVMCC_computeRhsJacobFastSingleSys("NumJacobFastSingleSys");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobFastSingleSys::FVMCC_ComputeRhsJacobFastSingleSys
(const std::string& name) :
  FVMCC_ComputeRhsJacobFastCoupling(name),
  _currLSSIdx(0)
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobFastSingleSys::~FVMCC_ComputeRhsJacobFastSingleSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobFastSingleSys::execute()
{
  CFLogDebugMin( "FVMCC_ComputeRhsJacobFastSingleSys::execute()" << "\n");

  //reset rhs to 0
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs = 0.0;

  const EquationSubSysDescriptor& eqSS =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  _currLSSIdx = eqSS.getEqSS();

  // set the equation subsystem descriptor
  const vector<CFuint>& currEqs = *_equations[_currLSSIdx];
  _nbEqs = currEqs.size();
  _start = currEqs[0];
  
  if (getMethodData().doComputeJacobian()) {
    _lss[_currLSSIdx]->getMatrix()->resetToZeroEntries();
  }
  
  // gradients and limiters are computed on all the variables at once
  _polyRec->computeGradients();
  _polyRec->computeLimiters();

  // also the nodal reconstruction from the cell centers to the vertices
  // is done at once
  if (_fluxData->hasDiffusiveTerm) {
    _nodalExtrapolator->extrapolateInAllNodes();
  }

  // set the list of faces
  vector<Common::SafePtr<TopologicalRegionSet> > trs =
    MeshDataStack::getActive()->getTrsList();

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

        // build the GeometricEntity
        geoData.idx = iFace;
        _currFace = geoBuilder->buildGE();

	// set the data for the FaceIntegrator
	setFaceIntegratorData();

	// extrapolate (and LIMIT, if the reconstruction is linear or more)
	// the solution in the quadrature points
	_polyRec->extrapolate(_currFace);

	// compute the physical data for each left and right reconstructed
	// state and in the left and right cell centers
	computeAndBackUpStatesData();

	const bool isBFace = _currFace->getState(1)->isGhost();
	if (!isBFace) {
	  _fluxSplitter->computeFlux(_flux);
	}
	else {
	  _currBC->computeFlux(_flux);
	}

	if (_fluxData->hasDiffusiveTerm) {
	  // reset to false the flag telling to freeze the diffusive
	  // coefficients
	  _diffVar->setFreezeCoeff(false);
	  // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
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
	    computeBoundaryJacobianTerm(_currLSSIdx);
	  }
	}
	
	// set off the perturbation flag
	_fluxData->isPerturb = false;
        geoBuilder->releaseGE();
      }
    }
  }

  if (getMethodData().isAxisymmetric() || getMethodData().hasSourceTerm()) {
    computeSourceTermContribution();
  }

  // compute dynamically the CFL
  getMethodData().getCFL()->setUpdateOn();
  getMethodData().getCFL()->update();
  getMethodData().getCFL()->setUpdateOff();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobFastSingleSys::computeJacobianTerm()
{
  if (_currFace->getState(0)->isParUpdatable() &&
      _currFace->getState(1)->isParUpdatable()) {
    FVMCC_ComputeRhsJacobFastCoupling::computeBothJacobTerms(_currLSSIdx);
  }
  else if (_currFace->getState(0)->isParUpdatable()) {
    FVMCC_ComputeRhsJacobFastCoupling::computeJacobTerm(_currLSSIdx, 0);
  }
  else if (_currFace->getState(1)->isParUpdatable()) {
    FVMCC_ComputeRhsJacobFastCoupling::computeJacobTerm(_currLSSIdx, 1);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobFastSingleSys::computeSourceTermContribution()
{
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  const bool isAxisymm = getMethodData().isAxisymmetric();
  const CFuint lastST = _stComputers->size() - 1;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  Common::SafePtr<GeometricEntityPool<CellTrsGeoBuilder> >
    geoBuilder = getMethodData().getCellTrsGeoBuilder();

  SafePtr<CellTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  CellTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  const vector<CFuint>& currEqs = *_equations[_currLSSIdx];

  // here it is supposed that if a source term is axisymmetric,
  // also all the others are
  // Moreover, this way of computing the source term fails if
  // the source term depends on more than one state
  for (CFuint ist = 0; ist < _stComputers->size(); ++ist) {
    // add a diagonal block for the source term
    for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
      // build the GeometricEntity
      geoData.idx = iCell;
      GeometricEntity *const currCell = geoBuilder->buildGE();

      State *const currState = currCell->getState(0);
      // consider only parallel updatable States
      if (currState->isParUpdatable()) {

	// compute the unperturbed source term (alll the components)
	(*_stComputers)[ist]->computeSource(currCell, _flux);

	// update the RHS
	if (isAxisymm) {
	  const RealVector& stateCoord = currState->getCoordinates();
	  for (CFuint iEq = 0; iEq < totalNbEqs; ++iEq) {
	    rhs(iCell, iEq, totalNbEqs) += getResFactor()*_flux[iEq];
	    // all the rhs has to be divided by r
	    // if there are more than 1 source term, this does
	    // this operation has to be done only ONCE after the
	    // last source term contribution has been added,
	    if (ist == lastST) {
	      rhs(iCell, iEq, totalNbEqs) /= std::abs(stateCoord[1]); // 1/r
	    }
	  }
	}
	else {
	  for (CFuint iEq = 0; iEq < totalNbEqs; ++iEq) {
	    rhs(iCell, iEq, totalNbEqs) += getResFactor()*_flux[iEq];
	  }
	}
	
	if (getMethodData().doComputeJacobian()) {
	  _bAcc[_currLSSIdx]->setRowColIndex(0, currState->getLocalID());
	  
	  // compute the source term corresponding to the perturbed state
	  for (CFuint iEq = 0; iEq < _nbEqs; ++iEq) {
	    const CFuint iVar = currEqs[iEq];
	    CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
	    // perturb the given component of the state vector
	    _numericalJacob->perturb(iVar, (*currState)[iVar]);
	    
	    // in the source term check the physical model to see if the
	    // number of equations and _start variable ID are the correct ones
	    (*_stComputers)[ist]->computeSource(currCell, _pertFlux);
	    
	    // compute the finite difference derivative of the source term
	    _numericalJacob->computeDerivative(_flux.slice(_start),
					       _pertFlux.slice(_start),
					       _fluxDiff.slice(_start));
	    
	    if (isAxisymm) {
	      const RealVector& stateCoord = currState->getCoordinates();
	      cf_assert(stateCoord[YY] > 0.0);
	      _fluxDiff.slice(_start) *= -getResFactor()/stateCoord[1];
	    }
	    else {
	      _fluxDiff.slice(_start) *= -getResFactor();
	    }
	    
	    _bAcc[_currLSSIdx]->addValues(0, 0, iEq, &_fluxDiff[_start]);
	    
	    // restore the unperturbed value
	    _numericalJacob->restore((*currState)[iVar]);
	  }
	  
	  // add the values in the jacobian matrix
	  _lss[_currLSSIdx]->getMatrix()->addValues(*_bAcc[_currLSSIdx]);
	  
	  // reset to zero the entries in the block accumulator
	  _bAcc[_currLSSIdx]->reset();
	}
      }

      // release the cell
      geoBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
