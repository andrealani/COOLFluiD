#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRHS.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "MathTools/MatrixInverter.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "FiniteVolume/DerivativeComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_ComputeRHS, CellCenterFVMData, FiniteVolumeModule>
FVMCC_computeRHSProvider("FVMCC");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRHS::FVMCC_ComputeRHS(const std::string& name) :
  CellCenterFVMCom(name),
  socket_cellFlag("cellFlag"),
  socket_normals("normals"),
  socket_faceAreas("faceAreas"),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_states("states"),
  socket_volumes("volumes"),
  socket_isOutward("isOutward"),
  socket_nstates("nstates"),
  socket_limiter("limiter"),
  socket_gstates("gstates"),
  socket_nodes("nodes"),
  _fluxSplitter(CFNULL),
  _diffusiveFlux(CFNULL),
  _reconstrVar(CFNULL),
  _diffVar(CFNULL),
  _polyRec(CFNULL),
  _nodalExtrapolator(CFNULL),
  _solutionToUpdateMatTrans(CFNULL),
  _updateToSolutionVecTrans(CFNULL),
  _stComputers(CFNULL),
  _eqFilters(CFNULL),
  _currFace(CFNULL),
  _currBC(CFNULL),
  _hasDiffusiveTerm(false),
  _isDiffusionActive(true),
  _faceIdx(0),
  _rMid(1.0),
  _upFactor(0.,2),
  _upStFactor(0.,2),
  _invr(1., 2),
  _flux(),
  _dFlux(),
  _rFlux(),
  _jacobDummy(),
  _invJacobDummy(),
  _source(),
  _sourceJacobian(),
  _stNumJacobIDs(),
  _stAnJacobIDs(),
  _sourceJacobOnCell(2),
  _fluxData(CFNULL),
  _tempUnitNormal(),
  _rExtraVars(),
  _inverter(CFNULL)  
{
  addConfigOptionsTo(this);

  _freezeDiffCoeff = false;
  setParameter("FreezeDiffCoeff",&_freezeDiffCoeff);

  _analyticalDiffJacob = false;
  setParameter("AnalyticDiffJacob",&_analyticalDiffJacob);
  
  _extrapolateInNodes = false;
  setParameter("FullNodalExtrapolation",&_extrapolateInNodes);
  
  _useAnalyticalMatrix = true;
  setParameter("useAnalyticalMatrix",&_useAnalyticalMatrix);
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRHS::~FVMCC_ComputeRHS()
{
  
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS::unsetup()
{
  // this could just be avoided: supposing that the RealVector* is set from an existing data
  for (CFuint i = 0; i < _rExtraVars.size(); ++i) {
    deletePtr(_rExtraVars[i]);
  }
  
  CellCenterFVMCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >
    ("FreezeDiffCoeff", "Flag forcing to freeze diffusive coefficients");

  options.addConfigOption< bool >
    ("AnalyticDiffJacob", "Compute analytically the jacobian of the diffusive fluxes");  
  
  options.addConfigOption< bool >
    ("FullNodalExtrapolation", "Extrapolate the solution to nodes consistently (normally it is not necessary)");

  options.addConfigOption< bool >
    ("useAnalyticalMatrix", "Flag telling if to use analytical matrix."); 
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS::configure ( Config::ConfigArgs& args )
{
  CellCenterFVMCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS::execute()
{
  CFTRACEBEGIN;
 
  CFLog(VERBOSE, "FVMCC_ComputeRHS::execute() START\n");
  
  initializeComputationRHS();
  
  // set the list of faces
  vector<SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trs.size();

  _faceIdx = 0;
  
  // no variable perturbation is needed in explicit residual computation
  getMethodData().setIsPerturb(false);
  
  // prepare the building of the faces
  Common::SafePtr<GeometricEntityPool<FaceCellTrsGeoBuilder> > geoBuilder = getMethodData().getFaceCellTrsGeoBuilder();
  geoBuilder->getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  FaceCellTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  vector<bool> zeroGrad(PhysicalModelStack::getActive()->getNbEq(), false);
  const bool hasSourceTerm = (getMethodData().isAxisymmetric() || getMethodData().hasSourceTerm());
  
  // this could be set during set up with no guarantee that it will be effective:
  // a MethodStrategy could set it to a different value afterwards, before entering here
  geoData.allCells = getMethodData().getBuildAllCells();
  
  const vector<string>& noBCTRS = getMethodData().getTRSsWithNoBC();
  SafePtr<CFMap<CFuint, FVMCC_BC*> > bcMap = getMethodData().getMapBC();
  
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];
    
    CFLog(VERBOSE, "TRS name = " << currTrs->getName() << "\n");
    // the faces on the boundary of the partition don't have to
    // be processed (their fluxes could give NaN)
    if (currTrs->getName() != "PartitionFaces" && currTrs->getName() != "InnerCells" && 
	!binary_search(noBCTRS.begin(), noBCTRS.end(), currTrs->getName())) {
      
      if (currTrs->hasTag("writable")) {
	_currBC = bcMap->find(iTRS);
	
	// set the flag telling if the ghost states have to be placed on the face itself
	_currBC->setPutGhostsOnFace();
	
        CFLog(VERBOSE, "BC name = " << _currBC->getName() << "\n");
	
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
      geoData.faces = currTrs;
      
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace, ++_faceIdx) {
        CFLogDebugMed( "iFace = " << iFace << "\n");
	
    	// reset the equation subsystem descriptor
	PhysicalModelStack::getActive()->resetEquationSubSysDescriptor();
	
	// build the GeometricEntity
        geoData.idx = iFace;
        _currFace = geoBuilder->buildGE();
	
	if (_currFace->getState(0)->isParUpdatable() || 
	    (!_currFace->getState(1)->isGhost() && _currFace->getState(1)->isParUpdatable())) {
	  
	  // set the data for the FaceIntegrator
	  setFaceIntegratorData();
	  
	  // cout << "states = " << _currFace->getState(0)->getLocalID()  << ", " <<  _currFace->getState(1)->getLocalID() << endl;
	  
	  // extrapolate (and LIMIT, if the reconstruction is linear or more)
	  // the solution in the quadrature points
	  _polyRec->extrapolate(_currFace);
	
	  // compute the physical data for each left and right reconstructed
	  // state and in the left and right cell centers
	  CFLog(DEBUG_MIN, "FVMCC_ComputeRHS::execute() => before computePhysicalData()\n");
	  computePhysicalData();
	  CFLog(DEBUG_MIN, "FVMCC_ComputeRHS::execute() => after computePhysicalData()\n");
	  
	  // a jacobian free method requires to re-compute the update coefficient every time the 
	  // residual is calculated to get F*v from the finite difference formula
	  // in particular the time dependent part of the residual depend on a updateCoeff
	  // that has to be up-to-date
	  getMethodData().setIsPerturb(false);
	  
	  // cout << "L = " << *_currFace->getState(0) << endl;
	  // cout << "R = " << *_currFace->getState(1) << endl;
	  
	  const bool isBFace = _currFace->getState(1)->isGhost();
	  
	  // this initialization is fundamental, especially for cases with coupling
	  // where some equation subsystems don't have convective terms
	  _flux = 0.; 
	  CFLog(DEBUG_MIN, "FVMCC_ComputeRHS::execute() => before conv computeFlux()\n");
	  if (!isBFace) {
	    _fluxSplitter->computeFlux(_flux);


	  }
	  else {
	    _currBC->computeFlux(_flux);
	  }
	  CFLog(DEBUG_MIN, "FVMCC_ComputeRHS::execute() => after conv computeFlux()\n");
	  // cout.precision(12);cout << currTrs->getName() << " C flux = " << _flux << endl;
	  
	  computeInterConvDiff();
	  
	  _isDiffusionActive = (*_eqFilters)[0]->filterOnGeo(_currFace);
	  
	  if (_hasDiffusiveTerm && _isDiffusionActive) {
	    // reset to false the flag telling to freeze the diffusive coefficients
	    _diffVar->setFreezeCoeff(false);
	    
	    // put virtual function here or parameter
	    if (_extrapolateInNodes) {
	      _nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
	    }
	    
	    // this initialization is fundamental, especially for cases with coupling
	    // where some equation subsystems don't have diffusive terms
	    _dFlux = 0.;
	    CFLog(DEBUG_MIN, "FVMCC_ComputeRHS::execute() => before diff computeFlux()\n");
	    _diffusiveFlux->computeFlux(_dFlux);
	    CFLog(DEBUG_MIN, "FVMCC_ComputeRHS::execute() => after diff computeFlux()\n");
	    _flux -= _dFlux;
	  }

	  /*if (_flux.size() == 9 && currTrs->getName() == "InnerFaces") {
	    cout.precision(12);cout << "SL[" << _currFace->getState(0)->getLocalID() << "] = [" << *_currFace->getState(0) << "]\n";
	    cout.precision(12);cout << "SR[" << _currFace->getState(1)->getLocalID() << "] = [" << *_currFace->getState(1) << "]\n";
	    cout.precision(12);cout << "["<< iFace << "] in " << currTrs->getName() << " C+D flux = " << _flux << endl; 
	    EXIT_AT(100);
	    }*/
	  
	  CFLogDebugMed("flux = " <<  _flux  << "\n");
	  
	  // compute the source term
	  if (hasSourceTerm) {
	    CFLog(DEBUG_MIN, "FVMCC_ComputeRHS::execute() => before computeSourceTerm()\n");
	    computeSourceTerm();
	    CFLog(DEBUG_MIN, "FVMCC_ComputeRHS::execute() => after computeSourceTerm()\n");
	  }
	  
	  // compute the contribution to the RHS
	  updateRHS();
	  // source term jacobians are only computed while processing internal faces 
	  CFLog(DEBUG_MIN, "FVMCC_ComputeRHS::execute() => before computeRHSJacobian()\n");
	  computeRHSJacobian();
	  CFLog(DEBUG_MIN, "FVMCC_ComputeRHS::execute() => after computeRHSJacobian()\n");
	}
	
	geoBuilder->releaseGE(); 
      }
    }
  }
    
  finalizeComputationRHS();
  
  /*const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  if (nbEqs == 9) {
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  ofstream fout("rhs.dat");
  for (CFuint iState = 0; iState < states.size(); ++iState) {
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      fout.precision(14); fout.setf(ios::scientific,ios::floatfield); fout << rhs(iState, iEq, nbEqs) << " ";
    }
    fout << endl;
  }
  }*/
  
  CFLog(VERBOSE, "FVMCC_ComputeRHS::execute() END\n");
  
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS::setup()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "FVMCC_ComputeRHS::setup() START\n");
  
  CellCenterFVMCom::setup();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _flux.resize(nbEqs,0.);
  _dFlux.resize(nbEqs,0.);
  _rFlux.resize(nbEqs,0.);
  _jacobDummy.resize(nbEqs, nbEqs, 0.);
  _invJacobDummy.resize(nbEqs, nbEqs, 0.);
  
  // set up the source terms
  _stComputers = getMethodData().getSourceTermComputer();
  const CFuint nbSourceTerms = _stComputers->size();
  _source.resize(2);
  _source[0].resize(nbSourceTerms);
  _source[1].resize(nbSourceTerms);
  for (CFuint k = 0; k < 2; ++k) {
    for (CFuint i = 0; i < nbSourceTerms; ++i) {
      _source[k][i].resize(nbEqs);
      _source[k][i] = 0.;
    }
  }
  
  _sourceJacobian.resize(2);
  _sourceJacobian[0].resize(nbSourceTerms);
  _sourceJacobian[1].resize(nbSourceTerms);
  for (CFuint k = 0; k < 2; ++k) {
    for (CFuint i = 0; i < nbSourceTerms; ++i) {
      _sourceJacobian[k][i].resize(nbEqs, nbEqs);
      _sourceJacobian[k][i] = 0.;
    }
  } 
  
  // differentiate between source terms with numerical or analytical jacobian 
  CFuint nbSourceTermsNumJacob = 0; 
  CFuint nbSourceTermsAnJacob = 0; 
  for (CFuint i = 0; i < nbSourceTerms; ++i) { 
    if ((*_stComputers)[i]->useAnalyticalJacob()) { 
      nbSourceTermsAnJacob++; 
    } 
    else { 
      nbSourceTermsNumJacob++; 
    } 
  } 
  
  // equation filters
  _eqFilters = getMethodData().getEquationFilters();
  cf_assert(_eqFilters->size() > 0);
  
  _stNumJacobIDs.reserve(nbSourceTermsNumJacob); 
  _stAnJacobIDs.reserve(nbSourceTermsAnJacob); 
  for (CFuint i = 0; i < nbSourceTerms; ++i) { 
    CFLog(DEBUG_MIN, "FVMCC_ComputeRHS::setup() => useAnalyticalJacob() for [" << 
	  i << "] => " << ((*_stComputers)[i]->useAnalyticalJacob()) << "\n"); 
    ((*_stComputers)[i]->useAnalyticalJacob()) ?  
      _stAnJacobIDs.push_back(i) : _stNumJacobIDs.push_back(i); 
  } 
  CFLog(DEBUG_MIN, "FVMCC_ComputeRHS::setup() => _stAnJacobIDs.size() = " << _stAnJacobIDs.size() << "\n");
  
  _tempUnitNormal.resize(PhysicalModelStack::getActive()->getDim(), 0.);
  _fluxSplitter = getMethodData().getFluxSplitter();
  
  _reconstrVar = (getMethodData().reconstructSolVars()) ?
    getMethodData().getSolutionVar() : getMethodData().getUpdateVar();
  
  _diffusiveFlux = getMethodData().getDiffusiveFluxComputer();
  _diffVar = getMethodData().getDiffusiveVar();
  _polyRec = getMethodData().getPolyReconstructor();
  _nodalExtrapolator = getMethodData().getNodalStatesExtrapolator();
  
  _solutionToUpdateMatTrans = getMethodData().getSolToUpdateInUpdateMatTrans();
  _updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  _solutionToUpdateMatTrans->setup(2);
  _updateToSolutionVecTrans->setup(2);
  
  getMethodData().getSolutionToLinearVecTrans()->setup(2); 
  getMethodData().getUpdateToReconstructionVecTrans()->setup(2);
  getMethodData().getReconstructionToUpdateVecTrans()->setup(2);
  getMethodData().getUpdateToSolutionInUpdateMatTrans()->setup(2);
  getMethodData().getJacobianLinearizer()->setMaxNbStates(2);
  
  // flag telling if a diffusive term has to be computed
  _hasDiffusiveTerm = !_diffusiveFlux->isNull();
  
  // set up the VarSets
  getMethodData().getUpdateVar()->setup();
  getMethodData().getSolutionVar()->setup();
  getMethodData().getDiffusiveVar()->setup();

  // change this!!!!
  _rExtraVars.resize(2);
  // this could just be avoided: supposing that the RelaVetor* is set from an existing data
  for (CFuint i = 0; i < 2; ++i) {
    _rExtraVars[i] = new RealVector(_reconstrVar->getExtraPhysicalVarsSize());
  }
  
  _inverter.reset(MatrixInverter::create(nbEqs, false));
  
  // set the builders of cells
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  
  // setup of the FaceCellTrsGeoBuilder
  SafePtr<FaceCellTrsGeoBuilder> geoBuilderPtr = getMethodData().getFaceCellTrsGeoBuilder()->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  geoBuilderPtr->setCellFlagSocket(socket_cellFlag);
  FaceCellTrsGeoBuilder::GeoData& geoData = getMethodData().getFaceCellTrsGeoBuilder()->getDataGE();
  geoData.cells = cells;
  
  // setup of the CellTrsGeoBuilder
  SafePtr<CellTrsGeoBuilder> cellBuilderPtr = getMethodData().getCellTrsGeoBuilder()->getGeoBuilder();
  cellBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  CellTrsGeoBuilder::GeoData& cellGeoData = getMethodData().getCellTrsGeoBuilder()->getDataGE();
  cellGeoData.trs = cells;
  
  CFLog(VERBOSE, "FVMCC_ComputeRHS::setup() END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS::updateRHS()
{
  if (getMethodData().isAxisymmetric()) {
    cf_assert(_currFace->nbNodes() == 2);
    const Node *const  node0 = _currFace->getNode(0);
    const Node *const  node1 = _currFace->getNode(1);
    
    // distance of the face mid point to the axis
    // _rMid = 0 on a centerline boundary face
    // _rMid == average y between the two nodes the face
    _rMid = 0.5*std::abs((*node0)[YY] + (*node1)[YY]);
    _invr[0] = 1./std::abs(_currFace->getState(0)->getCoordinates()[YY]);
    _invr[1] = 1./std::abs(_currFace->getState(1)->getCoordinates()[YY]);
  }
  
  const CFreal coeff = (getResFactor()*_rMid);
  _rFlux = coeff*_flux;
  
  // distribute the computed flux to the two neighbor states
  // of the corresponding face, subtracting the residual to the
  // current state and summing it to the neighbor state
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State *const firstState = _currFace->getState(0);
  const CFuint firstStateID = firstState->getLocalID();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    rhs(firstStateID, iEq, nbEqs) -= _rFlux[iEq]*_invr[0];
    
    CFLogDebugMed("rhs " << rhs(firstStateID, iEq, nbEqs) << "\n");
  }
  
  CFLogDebugMed("updateCoeff[firstStateID] = "
                << socket_updateCoeff.getDataHandle()[firstStateID] << "\n");
  
  // cout.precision(14); cout.setf(ios::scientific,ios::floatfield); cout << firstStateID << ", "   << socket_updateCoeff.getDataHandle()[firstStateID] << endl;
  
  State *const lastState = _currFace->getState(1);
  if (!lastState->isGhost()) {
    const CFuint lastStateID = lastState->getLocalID();
    cf_assert(firstStateID < lastStateID);
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      rhs(lastStateID, iEq, nbEqs) += _rFlux[iEq]*_invr[1];
      CFLogDebugMed("rhs " << rhs(lastStateID, iEq, nbEqs) << "\n");
    }
    
    CFLogDebugMed("updateCoeff[lastStateID] = "
	       << socket_updateCoeff.getDataHandle()[lastStateID]  << "\n");
    
    // cout.precision(14); cout.setf(ios::scientific,ios::floatfield); cout << lastStateID << ", "   << socket_updateCoeff.getDataHandle()[lastStateID] << endl;
  }
  
  // set the cell flags to true if the cells have already been used
  DataHandle<bool> cellFlag = socket_cellFlag.getDataHandle();
  for (CFuint i = 0; i < 2; ++i) { 
    if (!_currFace->getState(i)->isGhost()) {
      cellFlag[_currFace->getState(i)->getLocalID()] = true;
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS::setFaceIntegratorData()
{
  // normal points outward the first cell (== state) neighbor
  // set the current normal
  DataHandle< CFreal> normals = socket_normals.getDataHandle();
  cf_assert(_currFace != CFNULL);
  
  RealVector& unitNormal = getMethodData().getUnitNormal();
  const CFuint nbDim = unitNormal.size();
  const CFuint startID = _currFace->getID()*nbDim;
  const CFreal invArea = 1./socket_faceAreas.getDataHandle()[_currFace->getID()];
  for (CFuint i = 0; i < nbDim; ++i) {
    unitNormal[i] = normals[startID + i]*invArea;
  }
  
  getMethodData().getCurrentFace() = _currFace;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS::computeSourceTerm()
{
  CFTRACEBEGIN;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<bool> cellFlag = socket_cellFlag.getDataHandle();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // debugging -------------------------------------//
  
  for (CFuint iCell = 0; iCell < 2; ++iCell) { 

    if (_currFace->getState(iCell)->isGhost()){ 
      continue;
    }

    const CFuint cellID  = _currFace->getState(iCell)->getLocalID();
    _sourceJacobOnCell[iCell] = false;

    if (cellFlag[cellID] ){ 
      continue;
    }

    GeometricEntity *const currCell = _currFace->getNeighborGeo(iCell);
    CFreal invR = 1.0;
    if (getMethodData().isAxisymmetric()) {
      invR /= std::abs(currCell->getState(0)->getCoordinates()[YY]);
    }
    
    for (CFuint i = 0; i < _stAnJacobIDs.size(); ++i) { 
      const CFuint ist = _stAnJacobIDs[i]; 
      (*_stComputers)[ist]->setAnalyticalJacob(true);
      RealVector& source = _source[iCell][ist]; 
      source = 0.;
      (*_stComputers)[ist]->computeSource(currCell, source, _sourceJacobian[iCell][ist]); 
      
      CFLog(DEBUG_MED, "FVMCC_ComputeRHS::computeSourceTerm() => source = " << source << "\n"); 
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) { 
	rhs(cellID, iEq, nbEqs) += getResFactor()*source[iEq]*invR; 
      }
      
      (*_stComputers)[ist]->setAnalyticalJacob(false);
      _sourceJacobOnCell[iCell]= true;
    }
    
    for (CFuint i = 0; i < _stNumJacobIDs.size(); ++i) {  
      const CFuint ist = _stNumJacobIDs[i];
      RealVector& source = _source[iCell][ist];  
      source = 0.;
      (*_stComputers)[ist]->computeSource(currCell, source, _sourceJacobian[iCell][ist]);  
      
      CFLog(DEBUG_MED, "FVMCC_ComputeRHS::computeSourceTerm() => source = " << source << "\n"); 
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {  
	rhs(cellID, iEq, nbEqs) += getResFactor()*source[iEq]*invR;  
      }  
      
      cellFlag[cellID] = true;
      _sourceJacobOnCell[iCell]= true;
    }
  }
  
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS::transformResidual()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  RealVector tempRes(nbEqs);
  RealVector res(nbEqs);

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  
  for(CFuint iState = 0; iState < nbStates; ++iState) {
    // set and get the transformation matrix in the update variables
    const RealMatrix& matrix = (_useAnalyticalMatrix) ? 
      computeAnalyticalTransMatrix(*states[iState]) : computeNumericalTransMatrix(*states[iState]);
    
    // copy the rhs in a given temporary array
    const CFuint startID = iState*nbEqs;
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      tempRes[iEq] = rhs[startID + iEq];
    }
    
    // compute the transformed residual
    res = matrix*tempRes;
    
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      rhs[startID + iEq] = res[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& FVMCC_ComputeRHS::computeNumericalTransMatrix(State& state)
{
  // this first transformed state HAS TO BE stored,
  // since the returned pointer will change pointee after
  // the second call to transform()
  _rFlux = static_cast<RealVector&>(*_updateToSolutionVecTrans->transform(&state));
  
  NumericalJacobian& numJacob = getMethodData().getNumericalJacobian();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    // perturb the given component of the state vector
    numJacob.perturb(iVar, state[iVar]);
    
    const RealVector& tPertState = static_cast<RealVector&>
      (*_updateToSolutionVecTrans->transform(&state));
    
    // compute the finite difference derivative of the flux
    numJacob.computeDerivative(_rFlux, tPertState, _dFlux);
    _jacobDummy.setColumn(_dFlux,iVar);
    
    // restore the unperturbed value
    numJacob.restore(state[iVar]);
  }
  
  _inverter->invert(_jacobDummy, _invJacobDummy);
  return _invJacobDummy;
}

//////////////////////////////////////////////////////////////////////////////
      
vector<SafePtr<BaseDataSocketSink> > FVMCC_ComputeRHS::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_cellFlag);
  result.push_back(&socket_normals); 
  result.push_back(&socket_faceAreas); 
  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_states);
  result.push_back(&socket_volumes);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_nstates);
  result.push_back(&socket_limiter);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nodes);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS::initializeComputationRHS()
{
  // reset rhs to 0
  socket_rhs.getDataHandle() = 0.0;
  socket_cellFlag.getDataHandle() = false;
  
  for (CFuint i = 0; i < _eqFilters->size(); ++i) {
    (*_eqFilters)[i]->reset();
  }  
  
  // _polyRec->updateWeights();
  _polyRec->computeGradients();
  
  // extrapolate the solution from cell centers to all vertices
  _nodalExtrapolator->extrapolateInAllNodes();
}
      
//////////////////////////////////////////////////////////////////////////////
 
void FVMCC_ComputeRHS::computeRHSJacobian()
{
}
      
//////////////////////////////////////////////////////////////////////////////
   
void FVMCC_ComputeRHS::computeInterConvDiff()
{
}    

//////////////////////////////////////////////////////////////////////////////
      
void FVMCC_ComputeRHS::finalizeComputationRHS()
{
  if (getMethodData().isResidualTransformationNeeded()) {
    transformResidual();
  }
}

//////////////////////////////////////////////////////////////////////////////
 
 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
