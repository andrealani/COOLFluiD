#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCC_ComputeRhsJacobBlockDiag.hh"
#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_ComputeRhsJacobBlockDiag,
                      CellCenterFVMData,
                      FiniteVolumeModule>
fvmcc_ComputeRhsJacobBlockDiag("NumJacobBlockDiag");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobBlockDiag::FVMCC_ComputeRhsJacobBlockDiag(const std::string& name) :
  FVMCC_ComputeRHS(name),
  _lss(CFNULL),
  _numericalJacob(CFNULL),
  _pertFlux(),
  _fluxDiff(),
  _origState(),
  _acc(CFNULL),
  _bAcc(CFNULL),
  _cellNodes(),
  _axiFlux()
{
  addConfigOptionsTo(this);
// 	cout << "\n\n\n FVMCC_ComputeRhsJacobBlockDiag::FVMCC_ComputeRhsJacobBlockDiag \n\n\n";
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRhsJacobBlockDiag::~FVMCC_ComputeRhsJacobBlockDiag()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobBlockDiag::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobBlockDiag::execute()
{
  CFLogDebugMin( "FVMCC_ComputeRhsJacobBlockDiag::execute()" << "\n");
// 	cout <<"\n\n\n FVMCC_ComputeRhsJacobBlockDiag::execute() \n\n\n";

  //reset rhs to 0
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs = 0.0;

  // reset to zero all non zero entries in the jacobian
  if (getMethodData().doComputeJacobian()) {
    SafePtr<LSSMatrix> jacobMatrix = getMethodData().getLSSMatrix(0);
    jacobMatrix->resetToZeroEntries();
  }

  // add a computeVar !!!!
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
	computeStatesData();
	
	// a jacobian free method requires to re-compute the update coefficient every time the 
	// residual is calculated to get F*v from the finite difference formula
	// in particular the time dependent part of the residual depend on a updateCoeff
	// that has to be up-to-date
	_fluxData->isPerturb = false;
	
	const bool isBFace = _currFace->getState(1)->isGhost();
	if (!isBFace) {
	  _fluxSplitter->computeFlux(_flux);
	}
	else {
	  _currBC->computeFlux(_flux);
	}
		
	if (_firstOrderJacob) {
	  _fluxData->useFirstOrderJacob = _firstOrderJacob;
	  _fluxSplitter->computeFlux(_flux1ord);
	  _fluxData->useFirstOrderJacob = false;
	}
	
	if (_fluxData->hasDiffusiveTerm) {
	  // reset to false the flag telling to freeze the diffusive coefficients
	  _diffVar->setFreezeCoeff(false);
	  //_nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
	  _diffusiveFlux->computeFlux(_dFlux);
	  _flux -= _dFlux;
	  
	  if (_firstOrderJacob) {
            _flux1ord -= _dFlux;
          }
	}
	
	//compute the contribution to the RHS
	updateRHS();
	
	if (getMethodData().doComputeJacobian()) {
	  _fluxData->useFirstOrderJacob = _firstOrderJacob;
	  _diffVar->setFreezeCoeff(_freezeDiffCoeff);
	  
	  if (!isBFace) {
	    computeJacobianTerm();
	  }
	  else {
	    computeBoundaryJacobianTerm();
	  }
	}
	
	// set off the perturbation flag 
	_fluxData->isPerturb = false; 
	_fluxData->useFirstOrderJacob = false; 
	
	geoBuilder->releaseGE();
      }
    }
  }
  
  // reset the flag for freezing the transport properties
  _diffVar->setFreezeCoeff(false);

  if (getMethodData().isAxisymmetric() || getMethodData().hasSourceTerm()) {
    computeSourceTermContribution();
  }

  // compute dynamically the CFL
  if (getMethodData().doComputeJacobian()) {
    getMethodData().getCFL()->setUpdateOn();
    getMethodData().getCFL()->update();
    getMethodData().getCFL()->setUpdateOff();
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobBlockDiag::setup()
{
// 	cout <<"\n\n\n FVMCC_ComputeRhsJacobBlockDiag::setup \n\n\n" ;
  FVMCC_ComputeRHS::setup();

  _pertFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _fluxDiff.resize(PhysicalModelStack::getActive()->getNbEq());
  _origState.resize(PhysicalModelStack::getActive()->getNbEq());
  _axiFlux.resize(PhysicalModelStack::getActive()->getNbEq());

  // linear system solver
  _lss = getMethodData().getLinearSystemSolver()[0];

  // numerical jacobian
  _numericalJacob = &getMethodData().getNumericalJacobian();

  _acc.reset(_lss->createBlockAccumulator(2, 2, PhysicalModelStack::getActive()->getNbEq()));
  _bAcc.reset(_lss->createBlockAccumulator(1, 1, PhysicalModelStack::getActive()->getNbEq()));
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobBlockDiag::configure ( Config::ConfigArgs& args )
{
	//cout << "\n\n\n FVMCC_ComputeRhsJacobBlockDiag::configure BEGIN \n\n\n";
  FVMCC_ComputeRHS::configure(args);
	//cout << "\n\n\n FVMCC_ComputeRhsJacobBlockDiag::configure END \n\n\n";
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobBlockDiag::computeJacobianTerm()
{
	//cout <<"\n\n\n FVMCC_ComputeRhsJacobBlockDiag::computeJacobianTerm \n\n\n" ;
  if (_currFace->getState(0)->isParUpdatable() &&
      _currFace->getState(1)->isParUpdatable()) {
    computeBothJacobTerms();
  }
  else if (_currFace->getState(0)->isParUpdatable()) {
    computeJacobTerm(0);
  }
  else if (_currFace->getState(1)->isParUpdatable()) {
    computeJacobTerm(1);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobBlockDiag::computeBothJacobTerms()
{
	//cout << "\n\n\n FVMCC_ComputeRhsJacobBlockDiag::computeBothJacobTerms \n\n\n";
  // set on the perturbation flag
  _fluxData->isPerturb = true;
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
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
    _polyRec->extrapolate(_currFace, iVar, LEFT);

    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computeStatesData();

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
    _numericalJacob->computeDerivative(getJacobianFlux(),_pertFlux,_fluxDiff);
    // multiply for the residual factor
    _fluxDiff *= getResFactor();

    if (!getMethodData().isAxisymmetric()) {
      _acc->addValues(0, 0, iVar, &_fluxDiff[0]);

      // flux is opposite in sign for the other state
      _fluxDiff *= -1.;

      _acc->addValues(1, 0, iVar, &_fluxDiff[0]);
    }
    else {
      const CFreal invR0 = _rMid/(state0.getCoordinates())[YY];
      const CFreal invR1 = _rMid/(state1.getCoordinates())[YY];
      
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

    // linearization will be done in the flux splitter if needed
    _fluxSplitter->computeFlux(_pertFlux);

    if (_fluxData->hasDiffusiveTerm) {
      // _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);
      _pertFlux -= _dFlux;
    }

    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(getJacobianFlux(),
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
      const CFreal invR0 = _rMid/(state0.getCoordinates())[YY];
      const CFreal invR1 = _rMid/(state1.getCoordinates())[YY];
      
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
  getMethodData().getLSSMatrix(0)->addValues(*_acc);
  
  // reset to zero the entries in the block accumulator
  _acc->reset();
 }

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobBlockDiag::computeJacobTerm(const CFuint idx)
{
// 	cout << "\n\n\n FVMCC_ComputeRhsJacobBlockDiag::computeJacobTerm \n\n\n";
  cf_assert(_currFace->getState(idx)->isParUpdatable());
  
  // set on the perturbation flag
  _fluxData->isPerturb = true;
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
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
    _polyRec->extrapolate(_currFace, iVar, LEFT);
    
    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computeStatesData();
    
    // set the perturbed variable
    _fluxData->iPerturbVar = iVar;
    
    // linearization will be done in the flux splitter if needed
    _fluxSplitter->computeFlux(_pertFlux);
    
    if (_fluxData->hasDiffusiveTerm) {
      //_nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
      _diffusiveFlux->computeFlux(_dFlux);
      _pertFlux -= _dFlux;
    }
    
    // compute the finite difference derivative of the flux
    _numericalJacob->computeDerivative(getJacobianFlux(),_pertFlux,_fluxDiff);
    
    //multiply for the residual factor
    _fluxDiff *= factor;
    
    if (!getMethodData().isAxisymmetric() && state0.isParUpdatable()) { // state(idx) should be always parallel updatable
      _acc->addValues(idx, 0, iVar, &_fluxDiff[0]);
    }
    else if (state0.isParUpdatable()) { // state(idx) should be always parallel updatable
      const CFreal invR  = _rMid/(_currFace->getState(idx)->getCoordinates())[YY];
      _axiFlux = _fluxDiff*invR;
      _acc->addValues(idx, 0, iVar, &_axiFlux[0]);
    }
    
    // restore the unperturbed value
    _numericalJacob->restore(state0[iVar]);
    _polyRec->restoreValues(iVar, LEFT);
  }
  
  // second node contribution
  for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    CFLogDebugMax( "Perturbing iVar = " << iVar << "\n");
    
    // perturb the given component of the state vector
    _numericalJacob->perturb(iVar, state1[iVar]);
    // extrapolate (and LIMIT, if the reconstruction is linear or more)
    // the solution in the quadrature points
    _polyRec->extrapolate(_currFace, iVar, RIGHT);
    
    // compute the physical data for each left and right reconstructed
    // state and in the left and right cell centers
    computeStatesData();
    
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
    _numericalJacob->computeDerivative(getJacobianFlux(),
                                       _pertFlux,
                                       _fluxDiff);
    
    // multiply for the residual factor
    _fluxDiff *= factor;
    
    if (!getMethodData().isAxisymmetric() && state1.isParUpdatable()) { // state(idx) should be always parallel updatable
      _acc->addValues(idx, 1, iVar, &_fluxDiff[0]);
    }
    else if (state1.isParUpdatable()) { // state(idx) should be always parallel updatable
      const CFreal invR   = _rMid/(_currFace->getState(idx)->getCoordinates())[YY];
      _axiFlux = _fluxDiff*invR;
      _acc->addValues(idx, 1, iVar, &_axiFlux[0]);
    }
    
    // restore the unperturbed value
    _numericalJacob->restore(state1[iVar]);
    _polyRec->restoreValues(iVar, RIGHT);
  }
  
  // add the values in the jacobian matrix
  getMethodData().getLSSMatrix(0)->addValues(*_acc);
  
  // reset to zero the entries in the block accumulator
  _acc->reset();
}


//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobBlockDiag::computeBoundaryJacobianTerm()
{
  // set on the perturbation flag
  _fluxData->isPerturb = true;
  
  // set the index of the block corresponding to the current
  // state in the jacobian matrix
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State& currState = *_currFace->getState(0);
  State& ghostState = *_currFace->getState(1);
  cf_assert(ghostState.isGhost());
  
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
      
      if (_fluxData->hasDiffusiveTerm) {
	// _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
        _diffusiveFlux->computeFlux(_dFlux);
        _pertFlux -= _dFlux;
      }
      
      // compute the finite difference derivative of the flux
      _numericalJacob->computeDerivative(getJacobianFlux(),
                                         _pertFlux,
                                         _fluxDiff);

      // multiply for the residual factor
      _fluxDiff *= getResFactor();

      if (!getMethodData().isAxisymmetric()) {
        _bAcc->addValues(0, 0, iVar, &_fluxDiff[0]);
      }
      else {
	const CFreal invR0 = _rMid/(currState.getCoordinates())[YY];
	_axiFlux = _fluxDiff*invR0;
        _bAcc->addValues(0, 0, iVar, &_axiFlux[0]);
      }

      // restore the unperturbed value
      _numericalJacob->restore(currState[iVar]);

      // restore the original ghost state
      ghostState = _origState;
    }

    // add the values in the jacobian matrix
    getMethodData().getLSSMatrix(0)->addValues(*_bAcc);

    // reset to zero the entries in the block accumulator
    _bAcc->reset();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobBlockDiag::computeSourceTermContribution()
{
  /// @todo here it is supposed that if a source term is axisymmetric,
  ///       also all the others are
  if (getMethodData().doComputeJacobian()) {
    for (CFuint ist = 0; ist < _stComputers->size(); ++ist) {
     if (!(*_stComputers)[ist]->useAnalyticalJacob()) {
      computeSourceTermNumJacob(ist);
     }
     else {
       computeSourceTermAnalytJacob(ist);
     }
   }
 }
 else {
   FVMCC_ComputeRHS::computeSourceTerm();
 }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobBlockDiag::computeSourceTermAnalytJacob(CFuint ist)
{
  SafePtr<LSSMatrix> jacobMatrix = getMethodData().getLSSMatrix(0);

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const bool isAxisymm = getMethodData().isAxisymmetric();
  const CFuint lastST = _stComputers->size() - 1;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  CellTrsGeoBuilder::GeoData& geoData = _cellBuilder->getDataGE();

  // set on the analytical jacobian of the source term
  (*_stComputers)[ist]->setAnalyticalJacob(true);

  ///  @todo this way of computing the source term fails if
  ///        the source term depends on more than one state
  _flux = 0.0;
  
  //
  // static CFuint count = 0;
  // static std::vector<RealMatrix> st(socket_volumes.getDataHandle().size());
  //
    
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity *const currCell = _cellBuilder->buildGE();
    State *const currState = currCell->getState(0);
    
    if (currState->isParUpdatable()) {
      _bAcc->setRowColIndex(0, currState->getLocalID());
      
      // compute the unperturbed source term and its analytical jacobian
      (*_stComputers)[ist]->computeSource(currCell, _flux);
      
      // update the RHS
      if (isAxisymm) {
	const RealVector& stateCoord = currState->getCoordinates();
	const CFreal ovR = 1./stateCoord[YY]; // 1/r
	cf_assert(stateCoord[YY] > 0.0);
	
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  rhs(iCell, iEq, nbEqs) += getResFactor()*_flux[iEq];
	  // all the rhs has to be divided by r
	  // if there are more than 1 source term, this does
	  // this operation has to be done only ONCE after the
	  // last source term contribution has been added,
	  if (ist == lastST) {
	    rhs(iCell, iEq, nbEqs) *= ovR;
	  }
	}
      }
      else {
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  rhs(iCell, iEq, nbEqs) += getResFactor()*_flux[iEq];
	}
      }
      
      // get the jacobian matrix for the source term
      RealMatrix& sourceJacob = (*_stComputers)[ist]->getJacobMatrix();
      
      if (isAxisymm) {
	const RealVector& stateCoord = currState->getCoordinates();
	cf_assert(stateCoord[YY] > 0.0);
	sourceJacob *= -getResFactor()/stateCoord[YY];
      }
      else {
	sourceJacob *= -getResFactor();
      }
            
      _bAcc->addValues(0, 0, sourceJacob);
      
      // print the block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // _bAcc->print();

      // add the values in the jacobian matrix
      jacobMatrix->addValues(*_bAcc);

      // reset to zero the entries in the block accumulator
      _bAcc->reset();
    }

    // release the cell
    _cellBuilder->releaseGE();
  }

  // set off the analytical jacobian of the source term
  (*_stComputers)[ist]->setAnalyticalJacob(false);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRhsJacobBlockDiag::computeSourceTermNumJacob(CFuint ist)
{
  SafePtr<LSSMatrix> jacobMatrix = getMethodData().getLSSMatrix(0);

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const bool isAxisymm = getMethodData().isAxisymmetric();
  const CFuint lastST = _stComputers->size() - 1;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  CellTrsGeoBuilder::GeoData& geoData = _cellBuilder->getDataGE();

  ///  @todo this way of computing the source term fails if
  ///        the source term depends on more than one state
  // add a diagonal block for the source term

  _flux = 0.;
  _pertFlux = 0.0;

  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity *const currCell = _cellBuilder->buildGE();
    State *const currState = currCell->getState(0);
    
    if (currState->isParUpdatable()) {
      _bAcc->setRowColIndex(0, currState->getLocalID());
      // compute the unperturbed source term
      (*_stComputers)[ist]->computeSource(currCell, _flux);
      
      const vector<Node*>* const cellNodes = currCell->getNodes();
      const CFuint nbNodesInCell = cellNodes->size();
      
      // update the RHS
      if (isAxisymm) {
	const RealVector& stateCoord = currState->getCoordinates();
	const CFreal ovR = 1./stateCoord[YY]; // 1/r
	cf_assert(stateCoord[YY] > 0.0);
	
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  rhs(iCell, iEq, nbEqs) += getResFactor()*_flux[iEq];
	  // all the rhs has to be divided by r
	  // if there are more than 1 source term, this does
	  // this operation has to be done only ONCE after the
	  // last source term contribution has been added,
	  if (ist == lastST) {
	    rhs(iCell, iEq, nbEqs) *= ovR;
	  }
	}
      }
      else {
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  rhs(iCell, iEq, nbEqs) += getResFactor()*_flux[iEq];
	}
      }

      _fluxData->isPerturb = true;
      
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
	  _fluxDiff *= -getResFactor()/stateCoord[YY];
	}
	else {
	  _fluxDiff *= -getResFactor();
	}

	_bAcc->addValues(0, 0, iVar, &_fluxDiff[0]);

	// restore the unperturbed value
	_numericalJacob->restore((*currState)[iVar]);
      }

      // print the block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // _bAcc->print();

      // add the values in the jacobian matrix
      jacobMatrix->addValues(*_bAcc);
      
      // reset to zero the entries in the block accumulator
      _bAcc->reset();
      
      _fluxData->isPerturb = false;
    }
    
    // release the cell
    _cellBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD
