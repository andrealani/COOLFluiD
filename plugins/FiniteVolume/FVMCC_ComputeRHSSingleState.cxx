#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_ComputeRHSSingleState.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/VolumeIntegrator.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "FiniteVolume/ComputeDiffusiveFlux.hh"

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

MethodCommandProvider<FVMCC_ComputeRHSSingleState, CellCenterFVMData, FiniteVolumeModule>
FVMCC_computeRHSSingleStateProvider("FVMCCSingleState");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRHSSingleState::FVMCC_ComputeRHSSingleState(const std::string& name) :
  CellCenterFVMCom(name),
  socket_normals("normals"), 
  socket_faceAreas("faceAreas"),
  socket_states("states"),
  socket_isOutward("isOutward"),
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
  _hasDiffusiveTerm(false),
  _currFace(CFNULL),
  _cellBuilder(),
  _faceIdx(0),
  _rMid(0.0), 
  _flux(),
  _dFlux(),
  _rFlux(),
  _tempUnitNormal(), 
  _faceCoord(CFNULL),
  _rExtraVars()
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRHSSingleState::~FVMCC_ComputeRHSSingleState()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHSSingleState::unsetup()
{  
  deletePtr(_faceCoord);
  
  // this could just be avoided: supposing that the RelaVetor* is set from an existing data
  for (CFuint i = 0; i < _rExtraVars.size(); ++i) {
    deletePtr(_rExtraVars[i]);
  }
  
  CellCenterFVMCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHSSingleState::execute()
{
  CFTRACEBEGIN;
  
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  
  SafePtr <SpaceMethodData::PreconditionerData> pData = getMethodData().getPreconditionerData();
  
  const CFuint stateID = pData->currentStateID;
  
  const bool useBiggerIDs = pData->useBiggerStateIDs;
  const bool useAllIDs = pData->useAllStateIDs;
  
  RealVector& result = pData->result;
  
  /// @TODO check if this doesn't create problems
  // _polyRec->computeGradients();
  // _polyRec->computeLimiters();
  // extrapolate from cell centers to cell vertixes
  // _nodalExtrapolator->extrapolateInAllNodes();
  
  // no variable perturbation is needed in explicit residual computation
  getMethodData().setIsPerturb(true);
  
  vector<bool> zeroGrad(PhysicalModelStack::getActive()->getNbEq(), false);
  _polyRec->setZeroGradient(&zeroGrad);
  
  CellTrsGeoBuilder::GeoData& geoData = _cellBuilder.getDataGE();
  
  // build the GeometricEntityfv
  geoData.idx = stateID;
  
  // partition faces have to be skipped !!!!
  
  GeometricEntity *const currCell = _cellBuilder.buildGE();
  const vector<GeometricEntity*>& faces = *currCell->getNeighborGeos();
  const CFuint nbFaces = faces.size();
  
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    _currFace = faces[iFace];
    
    // you cannot ask an ID to a ghost state
    if (!_currFace->getState(RIGHT)->isGhost() && !_currFace->getState(LEFT)->isGhost()) {
      const CFuint faceID = _currFace->getID();
      const CFuint neighborID = (_currFace->getState(LEFT)->getLocalID() != stateID) ?
	_currFace->getState(LEFT)->getLocalID() : _currFace->getState(RIGHT)->getLocalID();
			
      if(((useBiggerIDs) && (neighborID > stateID)) || ((!useBiggerIDs) && (neighborID < stateID)) || (useAllIDs)) {
	//cout << neighborID << " ";
				
	_faceIdx = faceID;
				
	setFaceIntegratorData();
	
	if (static_cast<CFuint>(isOutward[faceID]) != stateID) {
	  getMethodData().getUnitNormal() *= -1.;
	}
	
	// extrapolate (and LIMIT, if the reconstruction is linear or more)
	// the solution in the quadrature points
	_polyRec->extrapolate(_currFace);
	
	
	// compute the physical data for each left and right reconstructed
	// state and in the left and right cell centeris
	computeStatesData();
	
	_fluxSplitter->computeFlux(_flux);
	
	if (_hasDiffusiveTerm) {
	  _diffusiveFlux->computeFlux(_dFlux);
	  _flux -= _dFlux;
	}
	
	// distance of the face mid point to the axis
	_rMid = 1.0; // _rMid = 0 on a centerline boundary face
	if (getMethodData().isAxisymmetric()) {
	  cf_assert(_currFace->nbNodes() == 2);
	  const Node *const  node0 = _currFace->getNode(0);
	  const Node *const  node1 = _currFace->getNode(1);
	  // _rMid == average y between the two nodes the face
	  _rMid *= 0.5*std::abs((*node0)[YY] + (*node1)[YY]);
	}
	
	result += (getResFactor()*_rMid)*_flux;	
	
	// 	cout << "result = " << result << endl;	
	// 	cout << "_flux = " << _flux << endl;
	// 	abort();
	
	//cout << "[" << stateID << "][" << neighborID << "] = "<< result << endl;
				
	CFLogDebugMed("flux = " <<  _flux  << "\n");
      }
    }
  }
  
  _cellBuilder.releaseGE();
	
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHSSingleState::setup()
{
  CFAUTOTRACE;
  
  CellCenterFVMCom::setup();
  
  _flux.resize(PhysicalModelStack::getActive()->getNbEq());
  _dFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _rFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempUnitNormal.resize(PhysicalModelStack::getActive()->getDim());
  _fluxSplitter = getMethodData().getFluxSplitter();
  
  // flag telling if a diffusive term has to be computed
  _hasDiffusiveTerm = !_diffusiveFlux->isNull();
  
  _reconstrVar = (getMethodData().reconstructSolVars()) ?
    getMethodData().getSolutionVar() : getMethodData().getUpdateVar();
  
  _diffusiveFlux = getMethodData().getDiffusiveFluxComputer();
  _diffVar = getMethodData().getDiffusiveVar();
  _polyRec = getMethodData().getPolyReconstructor();
  _nodalExtrapolator = getMethodData().getNodalStatesExtrapolator();
  
  _solutionToUpdateMatTrans = getMethodData().getSolToUpdateInUpdateMatTrans();
  _updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  _faceCoord = new Node(true);
    
  // change this!!!!
  _rExtraVars.resize(2);
  // this could just be avoided: supposing that the RelaVetor* is set from an existing data
  for (CFuint i = 0; i < 2; ++i) {
    _rExtraVars[i] = new RealVector(_reconstrVar->getExtraPhysicalVarsSize());
  }
  
  // set the builders of cells
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  
  _cellBuilder.setup();
  SafePtr<CellTrsGeoBuilder> geoBuilderPtr = _cellBuilder.getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  CellTrsGeoBuilder::GeoData& geoData = _cellBuilder.getDataGE();
  geoData.trs = cells;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHSSingleState::setFaceIntegratorData()
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

vector<SafePtr<BaseDataSocketSink> > FVMCC_ComputeRHSSingleState::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_normals);
  result.push_back(&socket_faceAreas);
  result.push_back(&socket_states);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nodes);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics
  
} // namespace COOLFluiD
