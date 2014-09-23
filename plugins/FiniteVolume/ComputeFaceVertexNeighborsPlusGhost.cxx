#include "ComputeFaceVertexNeighborsPlusGhost.hh"
#include "Framework/MeshData.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/CFMultiMap.hh"
#include "Framework/State.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeFaceVertexNeighborsPlusGhost,
			    ComputeStencil,
			    FiniteVolumeModule,
			    1>
computeFaceVertexNeighborsPlusGhostProvider("FaceVertexPlusGhost");

//////////////////////////////////////////////////////////////////////////////
      
//std::pair<CFuint,CFuint> ComputeFaceVertexNeighborsPlusGhost::GLOBAL_STATE_ID = std::pair<CFuint,CFuint>();

      
ComputeFaceVertexNeighborsPlusGhost::ComputeFaceVertexNeighborsPlusGhost(const std::string& name) :
   ComputeStencil(name)
 {
 }

//////////////////////////////////////////////////////////////////////////////

ComputeFaceVertexNeighborsPlusGhost::~ComputeFaceVertexNeighborsPlusGhost()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFaceVertexNeighborsPlusGhost::operator() ()
{  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  
  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();
  
  const CFuint nbElemTypes = elementType->size();
  
  // determine the size of the CFMultiMap to be created by counting
  // the number of nodes * number of elements per type
  CFuint nbPairsNodeState = 0;
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
    nbPairsNodeState += (*elementType)[iType].getNbNodes()*
      (*elementType)[iType].getNbElems();
  }
  
  // default behaviour: ghost states for all TRS are considered
  if (_trsNames.size() == 0) {
    vector< SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();
    const CFuint nbTRSs = trsList.size();
    for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) { 
      SafePtr<TopologicalRegionSet> currTrs = trsList[iTRS];
      if (currTrs->getName() != "PartitionFaces" && currTrs->getName() != "InnerCells") {
	if (currTrs->hasTag("writable")) {
	  _trsNames.push_back(currTrs->getName());
	}
      }
    }
  }
  
  for (CFuint iTRS = 0; iTRS < _trsNames.size(); ++iTRS) { 
    SafePtr<TopologicalRegionSet> bTrs = MeshDataStack::getActive()->getTrs(_trsNames[iTRS]);
    const CFuint nbTrsFaces = bTrs->getLocalNbGeoEnts();
    for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
      nbPairsNodeState += bTrs->getNbNodesInGeo(iFace);
    }
  }
    
  // allocate the memory for the map node-state
  typedef CFMultiMap<CFuint, State*> MapNodeState;
  MapNodeState mapNodeState(nbPairsNodeState);
  
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  // at first a multi-map with key = Node* and value = State*
  // referencing Node is created, by looping over all the nodes
  // inside the element-node connectivity
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    const CFuint nbNodesInCell = cells->getNbNodesInGeo(iCell);
    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
      mapNodeState.insert(cells->getNodeID(iCell, iNode), states[iCell]);
    }
  }
  
  GeometricEntityPool<FaceTrsGeoBuilder> geoBuilder;
  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder.getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  geoBuilder.setup();
  
  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();
  geoData.isBFace = true;
  
  // build the boundary faces
  for (CFuint iTRS = 0; iTRS < _trsNames.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> bTrs = MeshDataStack::getActive()->getTrs(_trsNames[iTRS]);
    cf_assert(bTrs->getName() != "PartitionFaces");  
    geoData.trs = bTrs;
    
    const CFuint nbTrsFaces = bTrs->getLocalNbGeoEnts();
    for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
      // build the GeometricEntity
      geoData.idx = iFace;
      GeometricEntity* const currFace = geoBuilder.buildGE();
      const CFuint nbFaceNodes = currFace->nbNodes();
      for (CFuint in = 0; in < nbFaceNodes; ++in) {
	mapNodeState.insert(currFace->getNode(in)->getLocalID(), currFace->getState(1));
      }
      
      geoBuilder.releaseGE();
    }  
  }
  
  mapNodeState.sortKeys(); // the sorting is fundamental !!
  assert(mapNodeState.getSize() == nbPairsNodeState);
  
  typedef MapNodeState::MapIterator mapIt;
  typedef vector<State*> SetOfNeighbors;
  SetOfNeighbors neighborList;
  neighborList.reserve(1000); // atomic number
  
  // AL: those are for debugging purposes
  static CFuint minStencil = 10000;
  static CFuint maxStencil = 0;
  static State* minState = 0;
  
  // Loop is made over all the cells (remember that iCell == iState !!!)
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    assert(iCell == states[iCell]->getLocalID());
    const CFuint nbNodesInCell = cells->getNbNodesInGeo(iCell);
    
    // loop over the nodes of the current cell
    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
      const CFuint nodeID = cells->getNodeID(iCell, iNode);
      // Get all the states referencing that node and store all of
      // them except the state inside the current cell (in fact,
      // you are looking for its neighbors)
      bool fo = false;
      pair<mapIt,mapIt> statesRefNode = mapNodeState.find(nodeID,fo);
      cf_assert(fo);
      
      for (mapIt stateInMap = statesRefNode.first;
	   stateInMap != statesRefNode.second;
	   ++stateInMap) {
	State *const state = stateInMap->second;
	
	// if the stateID of the neighbor is less than the current one,
	// the pair formed by this state and the neighbor has been
	// already inserted while processing the neighbor
	if (state != states[iCell]) {
	  //neighborList->insert(SetOfNeighbors::value_type(state));
	  
	  // here we need to keep the same exact insertion order
	  // easy solution suits this case where only few entries are required
          bool oldEntry = false;	
          for (CFuint i = 0; i < neighborList.size(); ++i) {
	    if (neighborList[i] == state) {
              oldEntry = true;
              break;
            }
	  }
          if (!oldEntry) neighborList.push_back(state);
	  // cout <<"state IN = " << state->getLocalID() << ", " << state->getCoordinates() << endl;
	}
      }
    }
    
    // remove duplicated entries
    
    
    // store the set of neighbor states in neighborStates
    const CFuint nbNeighStates = neighborList.size();
    stencil[iCell].resize(nbNeighStates);
    
    // cout << states[iCell]->getCoordinates() << " => " << nbNeighStates << endl;
    
    if (nbNeighStates < minStencil) {
      minStencil = min(minStencil, nbNeighStates);
      minState = states[iCell];
    }
    maxStencil = max(maxStencil, nbNeighStates);
    
    SetOfNeighbors::const_iterator it;
    CFuint in = 0;
    for (it = neighborList.begin(); it != neighborList.end(); ++it) {
      stencil[iCell][in] = *it;
      
      // cout <<"state OUT = " << (*it)->getLocalID() << ", " << (*it)->getCoordinates() << endl;
      in++;
    }
    
    neighborList.clear(); 
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume
    
  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
