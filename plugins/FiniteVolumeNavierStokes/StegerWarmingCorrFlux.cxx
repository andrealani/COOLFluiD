#include "FiniteVolumeNavierStokes/StegerWarmingCorrFlux.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/EulerTerm.hh"


#include "Framework/MapGeoToTrsAndIdx.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<StegerWarmingCorrFlux,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
stegerWarmingCorrFluxProvider("StegerWarmingCorr");

//////////////////////////////////////////////////////////////////////////////

StegerWarmingCorrFlux::StegerWarmingCorrFlux(const std::string& name) :
  StegerWarmingFlux(name), 
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nodes("nodes"),
  _sixth(1./6.),
  _third(1./3.),
  _k1(1. + _sixth - _third),
  _avStateL(CFNULL),
  _avStateR(CFNULL),
  _flagNormalFace(),
  _stencil(),
  _upStates(4)
{
  addConfigOptionsTo(this);

  _TWall = 0.;
  setParameter("TWall",&_TWall);
  
  _jacobDissip = 0.3;
  setParameter("jacobDissip",&_jacobDissip);
  
  _sigma = 0.5;
  setParameter("sigma",&_sigma); 
  
  _sigmaLim = 0.75;
  setParameter("sigmaLim",&_sigmaLim); 
  
  _maxNbNormalFaces = 10000;
  setParameter("maxNbNormalFaces",&_maxNbNormalFaces); 
  
  _useLimiter = false;
  setParameter("useLimiter",&_useLimiter);
  
  _useUpwindPolyRec = false;
  setParameter("useUpwindPolyRec",&_useUpwindPolyRec);
  
  _wallTrsNames = vector<std::string>();
  setParameter("wallTrsNames",&_wallTrsNames);
}

//////////////////////////////////////////////////////////////////////////////

StegerWarmingCorrFlux::~StegerWarmingCorrFlux()
{
  deletePtr(_avStateL);
  deletePtr(_avStateR);
}

//////////////////////////////////////////////////////////////////////////////

void StegerWarmingCorrFlux::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("TWall", "Wall temperature.");
  
  options.addConfigOption< CFreal >("sigma",
	"Coefficient controlling the pressure gradient-based sensor.");

  options.addConfigOption< CFreal >("sigmaLim",
	 "Coefficient controlling the pressure gradient-based limiter.");
 
  options.addConfigOption< CFreal >("jacobDissip",
        "Coefficient controlling the carbuncle fix.");
  
  options.addConfigOption< CFuint >("maxNbNormalFaces",
       "Maximum number of normal faces for which the carbuncle fix is inactive (if not specified will be set automatically).");
  
  options.addConfigOption< bool >("useLimiter",
      "Flag telling if the limiter must be used.");
   
  options.addConfigOption< bool >("useUpwindPolyRec",
       "Flag telling if the upwind reconstruction must be used.");
   
  options.addConfigOption< vector<std::string> >
     ("wallTrsNames", "Names of the wall TRSs.");
}	

//////////////////////////////////////////////////////////////////////////////

void StegerWarmingCorrFlux::compute(RealVector& result)
{
  CellCenterFVMData& data = this->getMethodData(); 
  GeometricEntity& face = *data.getCurrentFace();
  RealVector& pData = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->getPhysicalData();
  
  _updateVarSet->setExtraData(true);
  
  CFuint countSupersonic = 0;
  CFuint countSubsonic = 0;
  
  // compute physical data corresponding to the left state
  //  _updateVarSet->computePhysicalData(data.getCurrLeftState(), pData);
  _updateVarSet->computePhysicalData(*face.getState(LEFT), pData);
  const CFreal pL = pData[EulerTerm::P];
  const CFreal machL = pData[EulerTerm::V]/pData[EulerTerm::A];
  (machL >= 1.0) ? countSupersonic++ : countSubsonic++;
  
  // compute physical data corresponding to the right state
  //_updateVarSet->computePhysicalData(data.getCurrRightState(), pData);
  _updateVarSet->computePhysicalData(*face.getState(RIGHT), pData);
  const CFreal pR = pData[EulerTerm::P];
  const CFreal machR = pData[EulerTerm::V]/pData[EulerTerm::A];
  (machR >= 1.0) ? countSupersonic++ : countSubsonic++;
  
  // compute pressure based weights (Druguet, Candler, Nompelis)
  const CFreal gradP = _sigma*(pR - pL)/min(pR,pL);
  CFreal om = 0.5*(1./(gradP*gradP + 1.));
  CFreal oEminOm = 1. - om;
  
  CFreal eps = 0.0;
  bool isWallFace = false;
  
  if (_maxNbNormalFaces > 0) {
    eps = (_flagNormalFace[face.getID()]) ? 0.0 : _jacobDissip;
    
    // if you are at the wall, your face is in the list and it has a ghost right state
    isWallFace = (_flagNormalFace[face.getID()] && face.getState(RIGHT)->isGhost());
  }
  else {
    computeEps(pData, machL, machR, countSubsonic, countSupersonic, eps);
  }
  _solutionVarSet->setJacobDissipCoeff(eps);
  
  const State& left = *face.getState(LEFT);
  const State& right = *face.getState(RIGHT);
  
  // CHECK THIS
  if (isWallFace) {
    oEminOm = om = 0.5;
  }
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // compute A+ in the average left state
  *_avStateL = oEminOm*left + om*right;
    
  // ----------------- read temperature from file ---- //
  vector<Node*>& nodesInFace = *face.getNodes();
  const CFuint nbNodesInFace = nodesInFace.size();
  static Common::SafePtr<Framework::MapGeoToTrsAndIdx> m_mapGeoToTrs = MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
  
  SafePtr<TopologicalRegionSet> trs = m_mapGeoToTrs->getTrs(face.getID());
  bool isTrsWall = false;
  for (CFuint i =0; i < _wallTrsNames.size(); ++i) {
    if (trs->getName() ==  _wallTrsNames[i]) {
      isTrsWall = true;
      break;
    }
  }
   
  // AL: This allows for considering a given surface mapping of temperature 
  if (isTrsWall) {
    // build the mapTrs2NodalValues storage
    SafePtr<NodalStatesExtrapolator<CellCenterFVMData> > nse = 
      this->getMethodData().getNodalStatesExtrapolator();
    SafePtr<vector<NodalStatesExtrapolator<CellCenterFVMData>::MapTrs2NodalValues*> > 
      mapTrs2NodalValues = nse->getMapTrs2NodalValues();
    if (mapTrs2NodalValues->size() > 0) {
      RealVector& tWallArray = *(*mapTrs2NodalValues)[7]->find(&*trs); 
      CFMap<CFuint,CFuint>* mapNodeIDs = nse->getMapTrs2NodeIDs()->find(&*trs);
      // compute the temperature in the face centroid as the average of the nodal values
      _TWall = 0.0;
      for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
	const CFuint localNodeID = nodesInFace[iNode]->getLocalID();
	_TWall += tWallArray[mapNodeIDs->find(localNodeID)];
      }
      _TWall /= nbNodesInFace; 
    }
  }
  //---------------------
  
  if (isWallFace && _TWall > 0.0) {
    (*_avStateL)[nbEqs - 1] = (*_avStateL)[nbEqs - 2] = _TWall;
    for (CFuint i = 0; i < dim; ++i) {
      (*_avStateL)[nbEqs -3 -i] = 0.0;
    }
  }
  
  _updateVarSet->computePhysicalData(*_avStateL, pData);
  _solutionVarSet->splitJacobian(_jacobPlus,
				 _jacobDummy,
				 _eValues,
				 data.getUnitNormal());
  
  // compute A- in the average right state
  *_avStateR = om*left + oEminOm*right;
  
  if (isWallFace && _TWall > 0.0) {
    (*_avStateR)[nbEqs - 1] = (*_avStateR)[nbEqs - 2] = _TWall;
    for (CFuint i = 0; i < dim; ++i) {
      (*_avStateR)[nbEqs -3 -i] = 0.0;
    }
  }
  
  _updateVarSet->computePhysicalData(*_avStateR, pData);
  _solutionVarSet->splitJacobian(_jacobDummy,
				 _jacobMin,
				 _eValues,
				 data.getUnitNormal());
  
  _statesLR[0] = &data.getPolyReconstructor()->getCurrLeftState();
  _statesLR[1] = &data.getPolyReconstructor()->getCurrRightState();
  
  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  
  if (_useUpwindPolyRec) {
    const bool limitNeeded = (std::abs(pR - pL) > _sigmaLim*min(pR,pL));
    if (!_useLimiter || (_useLimiter && !limitNeeded)) {
      // _statesLR are constantly reconstructed
      upwindReconstruct();
    }
    else {  
      // _statesLR are constantly reconstructed
      _solutionStates = _updateToSolutionVarTrans->transformFromRefData(&pdata);
    }
  }
  else {
    _solutionStates = _updateToSolutionVarTrans->transformFromRefData(&pdata);
  }
  
  // compute Steger-Warming flux
  const State& stateL = *(*_solutionStates)[0];
  const State& stateR = *(*_solutionStates)[1];
  result = _jacobPlus*stateL + _jacobMin*stateR;
  
  // compute update coefficient
  if (!getMethodData().isPerturb()) {
    DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
    const CFreal faceArea = socket_faceAreas.getDataHandle()[face.getID()]/
      data.getPolyReconstructor()->nbQPoints();
    
    // left contribution to update coefficient
    CFreal maxEV =
      _updateVarSet->getMaxEigenValue(pdata[0], data.getUnitNormal());
    
    const CFuint leftID = face.getState(0)->getLocalID();
    updateCoeff[leftID] += max(maxEV, 0.)*faceArea;
    
    if (!face.getState(1)->isGhost()) {
      // right contribution to update coefficient
      
      _tempUnitNormal = -1.0*data.getUnitNormal();
      maxEV = _updateVarSet->getMaxEigenValue(pdata[1],_tempUnitNormal);
      const CFuint rightID = face.getState(1)->getLocalID();
      updateCoeff[rightID] += max(maxEV, 0.)*faceArea;
    }
  }
  
  _updateVarSet->setExtraData(false);
}

//////////////////////////////////////////////////////////////////////////////

void StegerWarmingCorrFlux::setup()
{
  StegerWarmingFlux::setup();
  
  _avStateL = new State();
  _avStateR = new State();
  
  if (_useUpwindPolyRec) {
    _updateToSolutionVarTrans->setup(4);
  }
  
  buildFaceBCData();
}

//////////////////////////////////////////////////////////////////////////////
      
void StegerWarmingCorrFlux::buildFaceBCData()
{ 
  setStencil();
  
  // this can only work serial!!
  if (_wallTrsNames.size() > 0) {
    _flagNormalFace.resize(this->socket_faceAreas.getDataHandle().size());
    _flagNormalFace.assign(_flagNormalFace.size(), false);
    
    // locally built geo builders. At this stage there could be 
    // something missing that doesn't let use the ones owned by the method data
    // it is safer to use local ones
    GeometricEntityPool<FaceTrsGeoBuilder> faceBuilder;
    faceBuilder.setup();
    faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
    FaceTrsGeoBuilder::GeoData& faceData = faceBuilder.getDataGE();
    
    GeometricEntityPool<CellTrsGeoBuilder> cellBuilder;
    cellBuilder.setup();
    cellBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
    CellTrsGeoBuilder::GeoData& cellData = cellBuilder.getDataGE();
    
    SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
      getTrs("InnerCells");
    CellTrsGeoBuilder::GeoData& geoData = cellBuilder.getDataGE();
    geoData.trs = cells;
    
    // ofstream fout("blFaces.dat");
    for (CFuint j = 0; j < _wallTrsNames.size(); ++j) {
      SafePtr<TopologicalRegionSet> wallTRS = 
	MeshDataStack::getActive()->getTrs(_wallTrsNames[j]);
      faceData.trs = wallTRS;
            
      CFuint cellID = 0;
      const CFuint nbTrsFaces = wallTRS->getLocalNbGeoEnts(); 
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {   
	faceData.isBFace = true;
	faceData.idx = iFace;
	GeometricEntity *const face = faceBuilder.buildGE();
	cellID = face->getState(0)->getLocalID();
	
	cf_assert(face->getState(1)->isGhost());
	CFuint faceID = face->getID();
	_flagNormalFace[faceID] = true;
	
	CFuint countFaces = 1;
	bool lastFace = false;
	while (!lastFace && (countFaces < _maxNbNormalFaces)) {
	  cellData.idx = cellID;
	  CFuint oppositeIFace = 0;
	  GeometricEntity *const cell = cellBuilder.buildGE();
	  const vector<GeometricEntity*>& cellFaces = *cell->getNeighborGeos();
	  const CFuint nbCellFaces =  cellFaces.size();
	  
	  for (CFuint f = 0; f < nbCellFaces; ++f) {
	    if (cellFaces[f]->getID() == faceID) {
	      oppositeIFace = getMethodData().getOppositeIFace
		(f, PhysicalModelStack::getActive()->getDim(), 
		 cell->nbNodes());
	      
	      countFaces++;
	      faceID = cellFaces[oppositeIFace]->getID();
	      _flagNormalFace[faceID] = true;
	      
	      const CFuint leftID = cellFaces[oppositeIFace]->getState(0)->getLocalID();
	      if (!cellFaces[oppositeIFace]->getState(1)->isGhost()) {
		const CFuint rightID = cellFaces[oppositeIFace]->getState(1)->getLocalID();
		cellID = (leftID == cellID) ? rightID : leftID;
	      }
	      else {
		// if you reach a face with a ghost state, you have reached the opposite boundary: stop!
		lastFace = true;
		
		// max number of faces is set to the total number of faces
		_maxNbNormalFaces = countFaces;
	      }

	      // fout << cell->getState(0)->getCoordinates() << endl;
	      break;
	    }
	  }
	  cellBuilder.releaseGE();
	}
	// fout << endl;
	faceBuilder.releaseGE();
	faceData.isBFace = false;
      } 
    }
    // fout.close();
  }
  
  cout << "StegerWarmingCorrFlux::buildFaceBCData() => _maxNbNormalFaces = " << _maxNbNormalFaces << endl;
}

//////////////////////////////////////////////////////////////////////////////

void StegerWarmingCorrFlux::setStencil()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  cf_assert(nbFaces > 0);
  _stencil.resize(nbFaces);
    
  GeometricEntityPool<FaceTrsGeoBuilder> faceBuilder;
  faceBuilder.setup();
  faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  FaceTrsGeoBuilder::GeoData& faceData = faceBuilder.getDataGE();
  
  GeometricEntityPool<CellTrsGeoBuilder> cellBuilder;
  cellBuilder.setup();
  cellBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  CellTrsGeoBuilder::GeoData& cellData = cellBuilder.getDataGE();
  
  SafePtr<TopologicalRegionSet> cells = 
    MeshDataStack::getActive()->getTrs("InnerCells");
  cellData.trs = cells;
  
  // set the list of faces
  vector<Common::SafePtr<TopologicalRegionSet> > trs = 
    MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trs.size();
  
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];
    faceData.trs = currTrs;
    
    // the faces on the boundary of the partition don't have to
    // be processed (their fluxes could give NaN)
    if (!currTrs->hasTag("partition") && currTrs->hasTag("face")) {
      if (currTrs->hasTag("writable")) {
	faceData.isBFace = true;
      }
      else {
        faceData.isBFace = false;
      }
      
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
        CFLogDebugMed( "iFace = " << iFace << "\n");
	
        // build the GeometricEntity
        faceData.idx = iFace;
        GeometricEntity *const face = faceBuilder.buildGE();
	const CFuint faceID = face->getID();
	_stencil[faceID].reserve((faceData.isBFace == false) ? 4 : 3);
	_stencil[faceID].push_back(face->getState(LEFT));
	_stencil[faceID].push_back(face->getState(RIGHT));
	
	// left cell
	cellData.idx = face->getState(LEFT)->getLocalID();
	CFuint oppositeIFace = 0;
	GeometricEntity* cell = cellBuilder.buildGE();
	const vector<GeometricEntity*>& cellFacesL = *cell->getNeighborGeos();
	const CFuint nbCellFacesL =  cellFacesL.size();
	for (CFuint f = 0; f < nbCellFacesL; ++f) {
	  if (cellFacesL[f]->getID() == faceID) {
	    oppositeIFace = getMethodData().getOppositeIFace
	      (f, PhysicalModelStack::getActive()->getDim(), 
	       cell->nbNodes());
	    
	    State *const leftState  = cellFacesL[oppositeIFace]->getState(0);
	    State *const rightState = cellFacesL[oppositeIFace]->getState(1);
	    _stencil[faceID].push_back((cell->getState(0) == leftState) ? 
				      rightState : leftState);
	    break;
	  }
	}
	cellBuilder.releaseGE();
	
	// if the current face is not a boundary face, look for the 
	// distance-1 neighbor of the face
	if (faceData.isBFace == false) {
	  cellData.idx = face->getState(RIGHT)->getLocalID();
	  cell = cellBuilder.buildGE();
	  
	  const vector<GeometricEntity*>& cellFacesR = *cell->getNeighborGeos();
	  const CFuint nbCellFacesR =  cellFacesR.size();
	  for (CFuint f = 0; f < nbCellFacesR; ++f) {
	    if (cellFacesR[f]->getID() == faceID) {
	      oppositeIFace = getMethodData().getOppositeIFace
		(f, PhysicalModelStack::getActive()->getDim(), cell->nbNodes());
	      
	      State *const leftState  = cellFacesR[oppositeIFace]->getState(0);
	      State *const rightState = cellFacesR[oppositeIFace]->getState(1);
	      _stencil[faceID].push_back((cell->getState(0) == leftState) ? 
					rightState : leftState);
	      break;
	    }
	  }
	  cellBuilder.releaseGE();
	}
	
	faceBuilder.releaseGE();
      }
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void StegerWarmingCorrFlux::computeEps(RealVector& pData,
				       const CFreal& machL,
				       const CFreal& machR, 
				       CFuint& countSubsonic,
				       CFuint& countSupersonic,
				       CFreal& eps)
{ 
  // compute the mach numbers in all the cells which share at least 
  // a vertex with the current face
  //  if (!data.isPerturb) {
  GeometricEntity& face = *this->getMethodData().getCurrentFace();
  State *const lstate = face.getState(0);
  State *const rstate = face.getState(1);
  vector<State*> vs;
  const CFuint nbNodesInFace = face.nbNodes();
  for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
    const CFuint nodeID = face.getNode(iNode)->getLocalID();
    const vector<State*>& s = getMethodData().getNodalStatesExtrapolator()->getNodalStateNeighbors(nodeID);
    
    for (CFuint is = 0; is < s.size(); ++is) {
      State *const currState = s[is];
      if (currState != lstate && currState != rstate) {
	bool isStateFound = false;
	for (CFuint i = 0; i < vs.size(); ++i) {
	  if (vs[i] == currState) {
	    isStateFound = true;
	    break;
	  }
	}
	
	if (!isStateFound) {
	  vs.push_back(currState);
	  _updateVarSet->computePhysicalData(*currState, pData);
	  const CFreal mach = pData[EulerTerm::V]/pData[EulerTerm::A];
	  (mach >= 1.0) ? countSupersonic++ : countSubsonic++;
	}
      }
    }
    
    // cf_assert((vs.size() > 0 && s.size() > 2) || (vs.size() == 0 && s.size() == 2));
    cf_assert(countSupersonic > 0 || countSubsonic > 0);
  }
  
  // set the jacobian dissipation coefficient to cure the carbuncle phenomenon
  // const CFreal eps2 = 0.3*(1.-2.*om);
  
  // eps = (countSupersonic > 0 && countSubsonic > 0) ? 0.3 : 0.3*(1.-2.*om);
  //    eps = _jacobDissip*(1.-2.*om);
  
  if (countSupersonic > 0 && countSubsonic > 0) {
    // eps == 0.3 if we are parallel to the shock
    eps = ((machL >= 1. && machR >= 1.) || (machL < 1. && machR < 1.)) ? 0.3 : 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StegerWarmingCorrFlux::upwindReconstruct()
{
  GeometricEntity& face = *getMethodData().getCurrentFace();
  cf_assert(face.getID() < _stencil.size());
  
  const vector<State*>& ss = _stencil[face.getID()]; 
  const bool isBFace = face.getState(RIGHT)->isGhost();
  _upStates.resize(ss.size());
  _upStates[0] = ss[0];
  _upStates[1] = ss[1];
  _upStates[2] = ss[2];
  
  if (!isBFace) {
    _upStates[3] = ss[3];
    cf_assert(_upStates.size() == 4);
    _solutionStates = _updateToSolutionVarTrans->transform(&_upStates);
    vector<State*>& cs = *_solutionStates;
    cf_assert(cs.size() == 4);
    *_avStateL = (*cs[0])*_k1 + (*cs[1])*_third - (*cs[2])*_sixth;
    *_avStateR = (*cs[0])*_third + (*cs[1])*_k1 - (*cs[3])*_sixth;
  }
  else { 
    cf_assert(_upStates.size() == 3);
    _solutionStates = _updateToSolutionVarTrans->transform(&_upStates);
    vector<State*>& cs = *_solutionStates;
    cf_assert(cs.size() == 4);
    *_avStateL = (*cs[0])*_k1 + (*cs[1])*_third - (*cs[2])*_sixth;
    *_avStateR = *_avStateL;
    // (*cs[0])*_third + (*cs[1])*(2.*_third);
    //here use weighted reconstruction 
    //(*cs[0])*d1 + (*cs[1])*d2;
  }

  cf_assert(ss.size() == _upStates.size());
  
  cf_assert(_solutionStates != CFNULL);  
  // override the first two solution states with the left and right values
  (*_solutionStates)[0]->copyData(static_cast<RealVector&>(*_avStateL));
  (*_solutionStates)[1]->copyData(static_cast<RealVector&>(*_avStateR));
}
      
//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StegerWarmingCorrFlux::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = StegerWarmingFlux::needsSockets();
  
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nodes);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
