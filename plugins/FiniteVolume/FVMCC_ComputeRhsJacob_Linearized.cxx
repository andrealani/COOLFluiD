#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRhsJacob_Linearized.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"
#include "Framework/VolumeIntegrator.hh"
#include "PolyReconstructorLin.hh"
#include "Framework/SubSystemStatus.hh"

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

MethodCommandProvider<FVMCC_ComputeRhsJacob_Linearized, CellCenterFVMData, FiniteVolumeModule> FVMCC_ComputeRhsJacob_Linearized("NumJacobLin");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacob_Linearized::FVMCC_ComputeRhsJacob_Linearized(const std::string& name) :
  FVMCC_ComputeRhsJacob(name),
  _fluxData2(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacob_Linearized::~FVMCC_ComputeRhsJacob_Linearized()
 {
 }

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob_Linearized::setFaceIntegratorData()
{
  // compute by shape function-based interpolation the
  // coordinates of the quadrature points

  getMethodData().getVolumeIntegrator()->
    computeCoordinatesAtQuadraturePointsOnGeoEnt(_currFace,_faceCoord);

  // normal points outward the first cell (== state) neighbor
  // set the current normal
  setFaceNormal();

  ///@todo the size of faceCoord should be not bigger than needed
  const CFuint nbQPoints = _nbQPointsInFace[_faceIdx];
  for (CFuint iq = 0; iq < nbQPoints; ++iq) {
    _fluxData->leftValues[iq]->getCoordinates()  = *_faceCoord[iq];
    _fluxData->rightValues[iq]->getCoordinates() = *_faceCoord[iq];
    _fluxData2->leftValues[iq]->getCoordinates()  = _fluxData->leftValues[iq]->getCoordinates();
    _fluxData2->rightValues[iq]->getCoordinates() = _fluxData->rightValues[iq]->getCoordinates();
  }

  // set the number of quadrature points in the data
  _fluxData->nbQPoints = nbQPoints;
  _fluxData->face = _currFace;
  _fluxData2->nbQPoints = nbQPoints;
  _fluxData2->face = _currFace;

}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob_Linearized::execute()
{
  CFLogDebugMin( "FVMCC_ComputeRhsJacob_Linearized::execute()" << "\n");

  //reset rhs to 0
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs = 0.0;

  SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();
  // reset to zero all non zero entries in the jacobian
  if(!getMethodData().isSysMatrixFrozen()) {
    jacobMatrix->resetToZeroEntries();
  }

  // add a computeVar !!!!
  _polyRecLin->computeGradients();
  _polyRecLin->computeLimiters();

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
	_polyRecLin->setZeroGradient(_currBC->getZeroGradientsFlags());
      }
      else {
        geoData.isBFace = false;
	_polyRecLin->setZeroGradient(&zeroGrad);
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
	_polyRecLin->extrapolate(_currFace);
	// compute the physical data for each left and right reconstructed
	// state and in the left and right cell centers
	computeStatesData2();

	const bool isBFace = _currFace->getState(1)->isGhost();
	if (!isBFace) {
	  _fluxSplitter->computeFlux(_flux);
	}
	else {
	  _currBC->computeFlux(_flux);
	}

	if (_fluxData->hasDiffusiveTerm) {
	  // reset to false the flag telling to freeze the
	  // diffusive coefficients
	  _diffVar->setFreezeCoeff(false);
	  // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
	  _diffusiveFlux->computeFlux(_dFlux);
	  _flux -= _dFlux;
	}

	//compute the contribution to the RHS
	updateRHS();

        // compute the jacobian matrix if it is not frozen
        if(!getMethodData().isSysMatrixFrozen()) {
          // set on the perturbation flag
          _fluxData->isPerturb = true;

	if (_freezeDiffCoeff)
	{
	  _diffVar->setFreezeCoeff(true);
	}

          if (!isBFace)
          {
            computeJacobianTerm();
          }
          else
          {
            computeBoundaryJacobianTerm();
          }

          // set off the perturbation flag
          _fluxData->isPerturb = false;
        }
        geoBuilder->releaseGE();
      }
    }
  }

  if (getMethodData().isAxisymmetric() || getMethodData().hasSourceTerm())
  {
    computeSourceTermContribution();
  }

  // compute dynamically the CFL
  getMethodData().getCFL()->setUpdateOn();
  getMethodData().getCFL()->update();
  getMethodData().getCFL()->setUpdateOff();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob_Linearized::setup()
{
  FVMCC_ComputeRhsJacob::setup();

  // Max number of quadrature points
  const CFuint maxNbQPts =
      getMethodData().getVolumeIntegrator()->getMaxIntegratorPattern().totalNbPts();

  // Max number of vertices
  const CFuint maxNbVertex = MeshDataStack::getActive()->Statistics().getMaxNbNodesInCell();

  // Resize fluxData vectors
  _fluxData2 = getMethodData().getFluxSplitterData2();

  // flag telling if a diffusive term has to be computed
  _fluxData2->hasDiffusiveTerm = !_diffusiveFlux->isNull();

  _fluxData2->normal.resize(PhysicalModelStack::getActive()->getDim());
  _fluxData2->unitNormal.resize(PhysicalModelStack::getActive()->getDim());
  _fluxData2->resizeGeoShapeFunctions(maxNbQPts, maxNbVertex);

  (_fluxData2->leftValues).resize(maxNbQPts);
  (_fluxData2->rightValues).resize(maxNbQPts);
  for (CFuint i = 0; i < maxNbQPts; ++i) {
    (_fluxData2->leftValues)[i] = new State();
    // create the node to store the position of the extrapolated state
    Node* node = new Node();
    node->setIsOnMesh(false);
    (_fluxData2->leftValues)[i]->setSpaceCoordinates(node);
  }

  for (CFuint i = 0; i < maxNbQPts; ++i) {
    (_fluxData2->rightValues)[i] = new State();
    // create the node to store the position of the extrapolated state
    Node* node = new Node();
    node->setIsOnMesh(false);
    (_fluxData2->rightValues)[i]->setSpaceCoordinates(node);
  }

  _fluxData2->faceNormals = &socket_normals;
  _fluxData2->isOutward   = &socket_isOutward;
  _fluxData2->faceAreas   = &_faceAreas;
  _fluxData2->nodalStates = &socket_nstates;
  _fluxData2->updateCoeff = &socket_updateCoeff;

  _rStates.resize(maxNbQPts*4);
  CFuint start = 0;
  for (CFuint i = 0; i < maxNbQPts; ++i) {
    _rStates[start]     = _fluxData->leftValues[i];
    _rStates[start + 1] = _fluxData->rightValues[i];
    start += 2;
  }
  for (CFuint i = start; i < 2*maxNbQPts; ++i) {
    _rStates[start]     = _fluxData2->leftValues[i];
    _rStates[start + 1] = _fluxData2->rightValues[i];
    start += 2;
  }


  // set the flux data
  _fluxSplitter->setFluxData2(_fluxData2);
  ///@todo add this!!
  ///_diffusiveFlux->setFluxData2(_fluxDataLinearized);
    
  _polyRecLin = getMethodData().getPolyReconstructor().d_castTo<PolyReconstructorLin>();
  _polyRecLin->setLeftValues2Ptr(&_fluxData2->leftValues);
  _polyRecLin->setRightValues2Ptr(&_fluxData2->rightValues);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob_Linearized::setFaceNormal()
{
  DataHandle< CFreal> normals = socket_normals.getDataHandle();
  cf_assert(_currFace != CFNULL);
  const CFuint nbDim = _fluxData->normal.size();
  const CFuint startID = _currFace->getID()*nbDim;
  const CFreal invArea = 1./_faceAreas[_currFace->getID()];
  for (CFuint i = 0; i < nbDim; ++i) {
    _fluxData->normal[i] = normals[startID + i];
    _fluxData->unitNormal[i] = _fluxData->normal[i]*invArea;
    _fluxData2->normal[i] = normals[startID + i];
    _fluxData2->unitNormal[i] = _fluxData2->normal[i]*invArea;
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob_Linearized::computeBothJacobTerms()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);

  _acc->setRowColIndex(0, state0.getLocalID());
  _acc->setRowColIndex(1, state1.getLocalID());

  const CFreal invR0 = _rMid/(state0.getCoordinates())[YY];
  const CFreal invR1 = _rMid/(state1.getCoordinates())[YY];

  // first node contribution
  for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");

    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, state0[iVar]);

    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRecLin->extrapolate(_currFace, iVar, LEFT);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computeStatesData2();

    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;

    // linearization will be done in the flux splitter if needed
    _fluxSplitter->computeFlux(_pertFlux);

    if (_fluxData->hasDiffusiveTerm) {
      // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);
      _pertFlux -= _dFlux;
    }

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
    _polyRecLin->restoreValues2(iVar, LEFT);

    // perturb the given component of the second state vector
    _numericalJacob->perturb(iVar, state1[iVar]);

    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRecLin->extrapolate(_currFace, iVar, RIGHT);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computeStatesData2();

    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;

    // linearization will be done in the flux splitter if needed
    _fluxSplitter->computeFlux(_pertFlux);

    if (_fluxData->hasDiffusiveTerm) {
      // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);
      _pertFlux -= _dFlux;
    }

    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(_flux,
                                       _pertFlux,
                                       _fluxDiff);

    // multiply for the residual factor
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
    _polyRecLin->restoreValues2(iVar, RIGHT);
  }

  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);

  // reset to zero the entries in the block accumulator
  _acc->reset();
 }

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob_Linearized::computeJacobTerm(const CFuint idx)
{
  cf_assert(_currFace->getState(idx)->isParUpdatable());

  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal invR   = _rMid/(_currFace->getState(idx)->getCoordinates())[YY];
  const CFreal factor = pow(-1.,static_cast<CFreal>(idx))*getResFactor();

  State& state0 = *_currFace->getState(0);
  State& state1 = *_currFace->getState(1);

  _acc->setRowColIndex(0, state0.getLocalID());
  _acc->setRowColIndex(1, state1.getLocalID());

  // first node contribution
  for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");

    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, state0[iVar]);

    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRecLin->extrapolate(_currFace, iVar, LEFT);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computeStatesData2();

    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;

    // linearization will be done in the flux splitter if needed
    _fluxSplitter->computeFlux(_pertFlux);

    if (_fluxData->hasDiffusiveTerm) {
      // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);
      _pertFlux -= _dFlux;
    }

    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(_flux,_pertFlux,_fluxDiff);

    //multiply for the residual factor
    _fluxDiff *= factor;

    if (!getMethodData().isAxisymmetric()) {
      _acc->addValues(idx, 0, iVar, &_fluxDiff[0]);
    }
    else {
      _axiFlux = _fluxDiff*invR;
      _acc->addValues(idx, 0, iVar, &_axiFlux[0]);
    }

    // restore the unperturbed value
    _numericalJacob->restore(state0[iVar]);
    _polyRecLin->restoreValues2(iVar, LEFT);
  }

  // second node contribution
  for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");

    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, state1[iVar]);
    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRecLin->extrapolate(_currFace, iVar, RIGHT);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computeStatesData2();

    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;

    // linearization will be done in the flux splitter if needed
    _fluxSplitter->computeFlux(_pertFlux);

    if (_fluxData->hasDiffusiveTerm) {
      // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);
      _pertFlux -= _dFlux;
    }

    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(_flux,
                                       _pertFlux,
                                       _fluxDiff);

    // multiply for the residual factor
    _fluxDiff *= factor;

    if (!getMethodData().isAxisymmetric()) {
      _acc->addValues(idx, 1, iVar, &_fluxDiff[0]);
    }
    else {
      _axiFlux = _fluxDiff*invR;
      _acc->addValues(idx, 1, iVar, &_axiFlux[0]);
    }

    // restore the unperturbed value
    _numericalJacob->restore(state1[iVar]);
    _polyRecLin->restoreValues2(iVar, RIGHT);
  }

  // add the values in the jacobian matrix
  _lss->getMatrix()->addValues(*_acc);

  // reset to zero the entries in the block accumulator
  _acc->reset();
}


//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacob_Linearized::computeBoundaryJacobianTerm()
{
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());

  const CFreal invR0 = _rMid/(currState.getCoordinates())[YY];

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
      _polyRecLin->extrapolate(_currFace);

      // compute the physical data for each left and right reconstructed
      // state and in the left and right cell centers
      computeStatesData2();

      // set the perturbed variable
      _fluxData->iPerturbVar = iVar;

      _currBC->computeFlux(_pertFlux);

      if (_fluxData->hasDiffusiveTerm) {
	// _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
        _diffusiveFlux->computeFlux(_dFlux);
        _pertFlux -= _dFlux;
      }

      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(_flux,
                                         _pertFlux,
                                         _fluxDiff);

      // multiply for the residual factor
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

void FVMCC_ComputeRhsJacob_Linearized::computeSourceTermContribution()
{
  SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
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

  /// @todo here it is supposed that if a source term is axisymmetric,
  ///       also all the others are
  for (CFuint ist = 0; ist < _stComputers->size(); ++ist) {
    ///  @todo this way of computing the source term fails if
    ///        the source term depends on more than one state
    // add a diagonal block for the source term
    for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
      // build the GeometricEntity
      geoData.idx = iCell;
      GeometricEntity *const currCell = geoBuilder->buildGE();

      State *const currState = currCell->getState(0);
      _bAcc->setRowColIndex(0, currState->getLocalID());
      // compute the unperturbed source term
      (*_stComputers)[ist]->computeSource(currCell, _flux);

      const vector<Node*>* const cellNodes = currCell->getNodes();
      const CFuint nbNodesInCell = cellNodes->size();

      if (currState->isParUpdatable()) {
        // update the RHS
        if (isAxisymm) {
          const RealVector& stateCoord = currState->getCoordinates();
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
            rhs(iCell, iEq, nbEqs) += getResFactor()*_flux[iEq];
            // all the rhs has to be divided by r
            // if there are more than 1 source term, this does
            // this operation has to be done only ONCE after the
            // last source term contribution has been added,
            if (ist == lastST) {
              rhs(iCell, iEq, nbEqs) /= std::abs(stateCoord[1]); // 1/r
            }
          }
        }
        else {
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
            rhs(iCell, iEq, nbEqs) += getResFactor()*_flux[iEq];
          }
        }

        // compute the source term corresponding to the perturbed state
        for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
          CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
          // perturb the given component of the state vector
          _numericalJacob->perturb(iVar, (*currState)[iVar]);

          // if it has diffusive flux, nodal reconstruction is done
          // starting from the perturbed state
          if (_fluxData->hasDiffusiveTerm) {
            _cellNodes.clear();
            // store the cells
            for (CFuint in = 0; in < nbNodesInCell; ++in) {
              _cellNodes.push_back((*cellNodes)[in]);
            }
            //_nodalExtrapolator->extrapolateInNodes(_cellNodes);
          }

          (*_stComputers)[ist]->computeSource(currCell, _pertFlux);
          // compute the finite difference derivative of the source term
          _numericalJacob->computeDerivative(_flux, _pertFlux, _fluxDiff);

          if (isAxisymm) {
            const RealVector& stateCoord = currState->getCoordinates();
            cf_assert(stateCoord[YY] > 0.0);
            _fluxDiff *= -getResFactor()/stateCoord[1];
          }
          else {
            _fluxDiff *= -getResFactor();
          }

          _bAcc->addValues(0, 0, iVar, &_fluxDiff[0]);

          // restore the unperturbed value
          _numericalJacob->restore((*currState)[iVar]);
        }

        // add the values in the jacobian matrix
        jacobMatrix->addValues(*_bAcc);

        // reset to zero the entries in the block accumulator
        _bAcc->reset();
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
