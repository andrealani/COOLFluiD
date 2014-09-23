#include "ComputeFaceEdgeNeighbors.hh"
#include "Framework/MeshData.hh"
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

Environment::ObjectProvider<ComputeFaceEdgeNeighbors,
			    ComputeStencil,
			    FiniteVolumeModule,
			    1>
computeFaceEdgeNeighborsProvider("FaceEdge");

//////////////////////////////////////////////////////////////////////////////

ComputeFaceEdgeNeighbors::ComputeFaceEdgeNeighbors(const std::string& name) :
   ComputeStencil(name)
 {
 }

//////////////////////////////////////////////////////////////////////////////

ComputeFaceEdgeNeighbors::~ComputeFaceEdgeNeighbors()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFaceEdgeNeighbors::operator() ()
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

  // allocate the memory for the map node-state
  typedef CFMultiMap<CFuint, State*> MapNodeState;
  MapNodeState mapNodeState(nbPairsNodeState);

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();

   // at first a multi-map with key = Node* and value = State*
  // referencing Node is created, by looping over all the nodes
  // inside the element-node connectivity
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    const CFuint nbNodesInCell = cells->getNbNodesInGeo(iCell);
    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
      mapNodeState.insert(cells->getNodeID(iCell, iNode),
			  states[iCell]);
    }
  }
  mapNodeState.sortKeys(); // the sorting is fundamental !!

  cf_assert(mapNodeState.getSize() == nbPairsNodeState);

  typedef MapNodeState::MapIterator mapIt;
  typedef set<State*, less<State*> > SetOfNeighbors;

  // Loop is made over all the cells (remember that iCell == iState !!!)
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    cf_assert(iCell == states[iCell]->getLocalID());

    const CFuint nbNodesInCell = cells->getNbNodesInGeo(iCell);
    SetOfNeighbors* neighborList = new SetOfNeighbors();
    // loop over the nodes of the current cell
    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
      const CFuint nodeID = cells->getNodeID(iCell, iNode);

      // Get all the states referencing that node and store all of
      // them except the state inside the current cell (in fact,
      // you are looking for its neighbors)
      bool fo = false;
      pair<mapIt,mapIt> statesRefNode = mapNodeState.find(nodeID, fo);
      cf_assert(fo);
      
      for (mapIt stateInMap = statesRefNode.first;
	   stateInMap != statesRefNode.second;
	   ++stateInMap) {
	State *const state = stateInMap->second;

	// if the stateID of the neighbor is less than the current one,
	// the pair formed by this state and the neighbor has been
	// already inserted while processing the neighbor
	const CFuint stateID = state->getLocalID();
	if (state != states[iCell] && (stateID > iCell)) {

	  // check if the considered state (cell) is an edge neighbor
	  // (== has at least 2 nodes in common with the current cell)
	  const CFuint nbNeighborNodes = cells->getNbNodesInGeo(stateID);
	  CFuint matchingNodes = 0;
	  for (CFuint in = 0; in < nbNodesInCell; ++in) {
	    // loop over the nodes of the neighbor cell
	    for (CFuint ik = 0; ik < nbNeighborNodes; ++ik) {
	      // add to the counter if match
	      // and break out of inner loop
	      if (cells->getNodeID(stateID,ik) == cells->getNodeID(iCell, in)) {
		++matchingNodes;
		break;
	      }
	    }
	  }

	  if (matchingNodes >= 2) {
	    neighborList->insert(SetOfNeighbors::value_type(state));
	  }
	}
      }
    }

    // store the set of neighbor states in neighborStates
    const CFuint nbNeighStates = neighborList->size();
    stencil[iCell].resize(nbNeighStates);

    SetOfNeighbors::iterator it;
    CFuint in = 0;
    for (it = neighborList->begin(); it != neighborList->end(); ++it) {
      stencil[iCell][in] = *it;
      in++;
    }
    deletePtr(neighborList);
  }
 }

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
