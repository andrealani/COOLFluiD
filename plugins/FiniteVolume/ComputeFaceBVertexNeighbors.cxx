#include "Common/CFMultiMap.hh"
#include "Framework/MeshData.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "FiniteVolume/ComputeFaceBVertexNeighbors.hh"
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

Environment::ObjectProvider<ComputeFaceBVertexNeighbors,
			    ComputeStencil,
			    FiniteVolumeModule,
			    1>
computeFaceBVertexNeighborsProvider("FaceBVertex");

//////////////////////////////////////////////////////////////////////////////

ComputeFaceBVertexNeighbors::ComputeFaceBVertexNeighbors(const std::string& name) :
   ComputeStencil(name)
 {
 }

//////////////////////////////////////////////////////////////////////////////

ComputeFaceBVertexNeighbors::~ComputeFaceBVertexNeighbors()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFaceBVertexNeighbors::operator() ()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();

  // create a sorted list of the boundary nodes
  vector<CFuint> bNodes;
  bNodes.reserve(nodes.size());

  // list of all the boundary nodes
  createBoundaryNodesList(bNodes);
  sort(bNodes.begin(), bNodes.end());

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  // list of the IDs of the cells that have at least a vertex on the boundary
  vector<CFuint> bCells;
  bCells.reserve(nbCells);

  // allocate the memory for the map nodeID - state pointer
  typedef CFMultiMap<CFuint, State*> MapNodeState;
  MapNodeState mapNodeState;

  GeometricEntityPool<CellTrsGeoBuilder> cellBuilder;
  SafePtr<CellTrsGeoBuilder> geoBuilderPtr = cellBuilder.getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  cellBuilder.setup();

  CellTrsGeoBuilder::GeoData& geoData = cellBuilder.getDataGE();
  geoData.trs = cells;

  // loop over ALL the cells to detected all the edges joining each
  // cell center and the centers of its FACE neighbors
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    State *const currState = states[iCell];
    const CFuint currStateID = currState->getLocalID();
    cf_assert(iCell == currStateID);

    // build the cell
    geoData.idx = iCell;
    GeometricEntity *const currCell = cellBuilder.buildGE();
    const GeomEntList* const neighborFaces = currCell->getNeighborGeos();
    const CFuint nbNeighborFaces = neighborFaces->size();

    bool isBCell = false;
    for (CFuint iFace = 0; iFace < nbNeighborFaces; ++iFace) {
      const GeometricEntity *const currFace =
	(*neighborFaces)[iFace];

      State *const neighborState =
	(currFace->getState(0) == currState) ?
	currFace->getState(1) : currFace->getState(0);

      if (!neighborState->isGhost()) {
        const CFuint neighborStateID = neighborState->getLocalID();
	cf_assert(neighborStateID < states.size());

        // if the stateID of the neighbor is less than the current one,
        // this means that the corresponding edge has already been
        // detected while processing the other stateID => skip it
        if (neighborStateID > currStateID) {
          stencil[currStateID].push_back(neighborState);
        }
      }
    }

    // distribute the state of the current cell to the nodes if
    // the node is a boundary node, then insert the current cell in the
    // list of boundary cells (== cells having at least ONE node on
    // the boundary)
    const vector<Node*> *const nodesInCell = currCell->getNodes();
    vector<Node*>::const_iterator itn;
    for (itn = nodesInCell->begin(); itn != nodesInCell->end(); ++itn) {
      if (binary_search(bNodes.begin(), bNodes.end(), (*itn)->getLocalID())) {
        mapNodeState.insert((*itn)->getLocalID(), currState);
        if (!isBCell) {
          isBCell = true;
          bCells.push_back(iCell);
        }
      }
    }

    // release the GeometricEntity
    cellBuilder.releaseGE();
  }

  mapNodeState.sortKeys(); // the sorting is fundamental !!

  typedef MapNodeState::MapIterator mapIt;
  //  typedef set<State*, less<State*> > SetOfNeighbors;

  // Loop is made only over the BOUNDARY cells
  vector<CFuint>::iterator itc;
  for (itc = bCells.begin(); itc != bCells.end(); ++itc) {
    const State *const currState = states[*itc];
    const CFuint currStateID = currState->getLocalID();
    cf_assert(*itc == currStateID);
    const CFuint nbNodesInCell = cells->getNbNodesInGeo(currStateID);

    for (CFuint in = 0; in < nbNodesInCell; ++in) {
      // get all the states referencing each node in the cell and
      // store the ones (!= current state) that are not yet
      // neighbors of the state
      const CFuint nodeID = cells->getNodeID(currStateID, in);
      if (binary_search(bNodes.begin(), bNodes.end(), nodeID)) {
	
	bool fo = false;
	pair<mapIt,mapIt> statesRefNode = mapNodeState.find(nodeID, fo);
	cf_assert(fo);
	
	for (mapIt stateInMap = statesRefNode.first;
	     stateInMap != statesRefNode.second;
	     ++stateInMap) {
	  State *const state = stateInMap->second;
	  if (state != currState) {
	    bool isAlreadyNeighbor = false;
	    std::vector<State*>::const_iterator its;
	    for (its = stencil[currStateID].begin();
		 its != stencil[currStateID].end(); ++its) {
	      if(*its == state) {
		isAlreadyNeighbor = true;
		break;
	      }
	    }
	    if (!isAlreadyNeighbor) {
	      stencil[currState->getLocalID()].push_back(state);
	    }
	  }
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFaceBVertexNeighbors::createBoundaryNodesList(vector<CFuint>& bNodes) const
{
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  // list of booleans to "colour" nodes and avoid multiple
  // insertions of the same nodes
  vector<bool> isInserted(nodes.size(), false);

  // loop over all the boundary TRSs to detect all the boundary nodes
  const vector< SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  vector< Common::SafePtr<TopologicalRegionSet> >::const_iterator itrs;
  for (itrs = trs.begin(); itrs != trs.end(); ++itrs) {
    const SafePtr<TopologicalRegionSet> currTrs = *itrs;
    if (currTrs->getName() != "InnerCells" &&
	currTrs->getName() != "InnerFaces") {

      // loop over all the faces
      const CFuint nbGeos = currTrs->getLocalNbGeoEnts();
      for (CFuint iGeo = 0; iGeo < nbGeos; ++iGeo) {
	const CFuint nbNodesInGeo = currTrs->getNbNodesInGeo(iGeo);
	for (CFuint iNode = 0; iNode < nbNodesInGeo; ++iNode) {
	  const CFuint nodeID = currTrs->getNodeID(iGeo, iNode);
	  if (!isInserted[nodeID]) {
	    bNodes.push_back(nodeID);
	    isInserted[nodeID] = true;
	  }
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
