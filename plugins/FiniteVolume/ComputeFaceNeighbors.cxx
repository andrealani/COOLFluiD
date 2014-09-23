#include "ComputeFaceNeighbors.hh"
#include "Framework/MeshData.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/CFMultiMap.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
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

Environment::ObjectProvider<ComputeFaceNeighbors,
			    ComputeStencil,
			    FiniteVolumeModule,
			    1>
computeFaceNeighborsProvider("Face");

//////////////////////////////////////////////////////////////////////////////

ComputeFaceNeighbors::ComputeFaceNeighbors(const std::string& name) :
   ComputeStencil(name)
 {
 }

//////////////////////////////////////////////////////////////////////////////

ComputeFaceNeighbors::~ComputeFaceNeighbors()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFaceNeighbors::operator() ()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();

  SafePtr<MapGeoToTrsAndIdx> mapGeoToTrs = MeshDataStack::getActive()->
    getMapGeoToTrs("MapFacesToTrs");

  SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  SafePtr<ConnectivityTable<CFuint> > cellFaces =
    MeshDataStack::getActive()->getConnectivity("cellFaces");

  const CFuint nbElems = cells->getLocalNbGeoEnts();

  // loop over ALL the cells to detected all the edges joining each
  // cell center and the centers of its FACE neighbors
  for (CFuint iElem = 0; iElem < nbElems; ++iElem) {
    State *const currState = states[iElem];
    const CFuint currStateID = currState->getLocalID();

    cf_assert(iElem == currStateID);

    const CFuint nbNeighborFaces = cellFaces->nbCols(iElem);
    for (CFuint iFace = 0; iFace < nbNeighborFaces; ++iFace) {
      const CFuint faceID = (*cellFaces)(iElem, iFace);
      const CFuint faceIdx = mapGeoToTrs->getIdxInTrs(faceID);
      SafePtr<TopologicalRegionSet> trs = mapGeoToTrs->getTrs(faceID);
      const bool isBFace = mapGeoToTrs->isBGeo(faceID);
      State *const state0 = states[trs->getStateID(faceIdx,0)];
      State* neighborState = CFNULL;
      const CFuint neighID = trs->getStateID(faceIdx,1);
      if (state0 == currState) {
	// first state of the face == current one
	// the second state for boundary faces is stored in the
	// gstates handle
	neighborState = (!isBFace) ?
	  states[neighID] : gstates[neighID];
      }
      else {
	neighborState = states[neighID];
      }

      // neglect the neighbor state if it is a ghost one
      if (!isBFace) {
	const CFuint neighborStateID = neighborState->getLocalID();
	// if the stateID of the neighbor is less than the current one,
	// this means that the corresponding edge has already been
	// detected while processing the other stateID => skip it
	if (neighborStateID > currStateID) {
	  stencil[currStateID].push_back(neighborState);
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
