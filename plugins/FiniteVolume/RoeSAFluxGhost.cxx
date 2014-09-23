#include "RoeSAFluxGhost.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/BaseTerm.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/CellTrsGeoBuilder.hh"

#include <fstream>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RoeSAFluxGhost,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
roeSAFluxGhostProvider("RoeSAGhost");
      
//////////////////////////////////////////////////////////////////////////////

RoeSAFluxGhost::RoeSAFluxGhost(const std::string& name) :
  RoeFlux(name),
  m_cellBuilder(),
  _isPartitionFace(),
  _dataLeftState(),
  _dataRightState(),
  _tmpEv(),
  _normal()
{
  addConfigOptionsTo(this);
  _entropyFixID = 0;
  setParameter("entropyFixID",&_entropyFixID);
}

//////////////////////////////////////////////////////////////////////////////

RoeSAFluxGhost::~RoeSAFluxGhost()
{
}

//////////////////////////////////////////////////////////////////////////////

void RoeSAFluxGhost::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("entropyFixID","ID of the entropy correction type.");
}

//////////////////////////////////////////////////////////////////////////////

void RoeSAFluxGhost::setAbsEigenValues()
{  
  //compute eigen values of the left state
  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  
  updateVarSet->computeEigenValues(pdata[0], 
				   unitNormal,
				   _leftEvalues);
  
  //compute eigen values of the right state
  updateVarSet->computeEigenValues(pdata[1],
				   unitNormal,
				   _rightEvalues);
  
  //compute eigen values of the left state
  GeometricEntity& face = *getMethodData().getCurrentFace();
  CFreal etaSA = 0.0;
  
  // 2D case
  // build the left cell
  if (!face.getState(1)->isGhost()) {
    updateEtaPA(face.getNeighborGeo(LEFT), _leftEvalues, etaSA);
    updateEtaPA(face.getNeighborGeo(RIGHT), _rightEvalues, etaSA);
  }
  else {
    // to be fixed
    updateEtaPAGhost(_leftEvalues, _rightEvalues, etaSA);
  }
  
  etaSA *= 0.5;
  assert(std::abs(etaSA) > 0.);
  
//   const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint i = 0; i < nbEqs; ++i) {
    const CFreal absEv = std::abs(_eValues[i]);
    
    switch (_entropyFixID) {
    case 1:
      _absEvalues[i] = (absEv >= 2.*etaSA) ? absEv : 0.25*absEv*absEv/etaSA + etaSA; 
      break;
    case 2:
      _absEvalues[i] = max(absEv, etaSA);
      break;
    case 3:
      _absEvalues[i] = absEv + etaSA;
      break;
    default:
      _absEvalues[i] = absEv + etaSA;
      break;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RoeSAFluxGhost::setup()
{
  RoeFlux::setup();
    
  m_cellBuilder.setup();
  m_cellBuilder.setDataSockets(socket_states, socket_gstates, socket_nodes);
  CellTrsGeoBuilder::GeoData& cellGeoData = m_cellBuilder.getDataGE();
  cellGeoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");
  
  PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->resizePhysicalData(_dataLeftState);
  PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->resizePhysicalData(_dataRightState);
  _tmpEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _normal.resize(PhysicalModelStack::getActive()->getDim());
  
  getMethodData().setBuildAllCells(true);  
  
  Common::SafePtr<TopologicalRegionSet> pFaces = MeshDataStack::getActive()->
    getTrs("PartitionFaces");
  const CFuint nbPFaces = pFaces->getLocalNbGeoEnts();
  _isPartitionFace.resize(MeshDataStack::getActive()->Statistics().getNbFaces(), false);
  for (CFuint i = 0; i < nbPFaces; ++i) {
    _isPartitionFace[pFaces->getLocalGeoID(i)] = true;
  }
  
}
      
//////////////////////////////////////////////////////////////////////////////

void RoeSAFluxGhost::updateEtaPA(GeometricEntity *const cell, 
			    const RealVector& eValues,
			    CFreal& etaSA) 
{
  GeometricEntity& currFace = *getMethodData().getCurrentFace();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  
  const GeomEntList *const faces = cell->getNeighborGeos();
  const CFuint nbFaces = faces->size();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    GeometricEntity *const neighborFace = (*faces)[iFace];
    // only faces sharing nodes with the orginal face are considered
    if (!_isPartitionFace[neighborFace->getID()]) {
      if (neighborFace->sharesNodesWith(currFace)) {
	State *const rs = (neighborFace->getState(0) == cell->getState(0)) ? 
	  neighborFace->getState(1) : neighborFace->getState(0);
	assert(rs != NULL);
	
	const CFuint startID = neighborFace->getID()*dim;
	assert(neighborFace->getID() < socket_faceAreas.getDataHandle().size());
	const CFreal invArea = 1./socket_faceAreas.getDataHandle()[neighborFace->getID()];
	
	for (CFuint i = 0; i < dim; ++i) {
	  _normal[i] = normals[startID + i]*invArea;  
	}	
	
	updateVarSet->computePhysicalData(*rs, _dataRightState);  
	updateVarSet->computeEigenValues(_dataRightState,_normal, _tmpEv);
	
	for (CFuint i = 0; i < nbEqs; ++i) {
	  etaSA = max(etaSA, std::abs(_tmpEv[i] - eValues[i]));       
	}
      }
    }
  }
}	
      
//////////////////////////////////////////////////////////////////////////////

void RoeSAFluxGhost::updateEtaPAGhost(const RealVector& lEvalues,
				      const RealVector& rEvalues,
				      CFreal& etaSA) 
{
  GeometricEntity& currFace = *getMethodData().getCurrentFace();
  GeometricEntity& cell = *currFace.getNeighborGeo(LEFT);
  State *const rState = currFace.getState(RIGHT);
  cf_assert(rState->isGhost());
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  
  CellTrsGeoBuilder::GeoData& geoData = m_cellBuilder.getDataGE();
  
  const GeomEntList *const faces = cell.getNeighborGeos();
  const CFuint nbFaces = faces->size();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    GeometricEntity *const neighborFace = (*faces)[iFace];
    // only faces sharing nodes with the original face are considered
    if (!_isPartitionFace[neighborFace->getID()]) {
      if (neighborFace->sharesNodesWith(currFace) && !neighborFace->getState(RIGHT)->isGhost()) {
	// here consider each of the two neighbors of the left cell
	State *const rs = (neighborFace->getState(0) == cell.getState(0)) ? 
	  neighborFace->getState(1) : neighborFace->getState(0);
	assert(rs != NULL);
	assert(!rs->isGhost());
	
	const CFuint startID = neighborFace->getID()*dim;
	assert(neighborFace->getID() < socket_faceAreas.getDataHandle().size());
	const CFreal invArea = 1./socket_faceAreas.getDataHandle()[neighborFace->getID()];
	
	for (CFuint i = 0; i < dim; ++i) {
	  _normal[i] = normals[startID + i]*invArea;  
	}	
	
	updateVarSet->computePhysicalData(*rs, _dataRightState);  
	updateVarSet->computeEigenValues(_dataRightState,_normal, _tmpEv);
	for (CFuint i = 0; i < nbEqs; ++i) {
	  etaSA = max(etaSA, std::abs(_tmpEv[i] - lEvalues[i]));       
	}
	
	// here create the cell corrisponding to the neighbor state of current left state 
	// and, while looping through its faces,  look for the only other face sharing
	// at least one node with the initial face: that's a ghost we have to consider 
	// in our stencil 
	const CFuint cellID = rs->getLocalID();
	geoData.idx = cellID;
	
	cout << "cellID = " << cellID << endl;
	GeometricEntity *const nCell = m_cellBuilder.buildGE();
	cout << "nCell = " << nCell << endl;
	
	const GeomEntList *const nfaces = nCell->getNeighborGeos();
	const CFuint nbnFaces = nfaces->size();
	for (CFuint f = 0; f < nbnFaces; ++f) {
	  GeometricEntity *const nf = (*nfaces)[f];
	  if (!_isPartitionFace[nf->getID()]) {
	    if (nf->sharesNodesWith(currFace) && nf->getID() != neighborFace->getID()) {
	      cf_assert(nf->getState(0)->isGhost() || nf->getState(1)->isGhost());
	      State *const rs2 = (nf->getState(0)->isGhost()) ?  nf->getState(0) : nf->getState(1);
	      assert(rs2->isGhost()); // careful: ghost states do not have valid IDs
	      
	      // here we define a fictitious normal between rs2 and rState, taking the direction joining the two points
	      _normal = rs2->getCoordinates() - rState->getCoordinates();
	      _normal.normalize();
	      
	      updateVarSet->computePhysicalData(*rs2, _dataRightState);  
	      updateVarSet->computeEigenValues(_dataRightState,_normal, _tmpEv);
	      for (CFuint i = 0; i < nbEqs; ++i) {
		etaSA = max(etaSA, std::abs(_tmpEv[i] - rEvalues[i]));       
	      }
	    }
	  }
	}
	
	m_cellBuilder.releaseGE();
      }
    }
  }
}	
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
