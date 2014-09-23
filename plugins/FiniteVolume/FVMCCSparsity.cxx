#include "Environment/ObjectProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/State.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCCSparsity.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FVMCCSparsity,
               GlobalJacobianSparsity,
               FiniteVolumeModule>
aFVMCCSparsityProvider("FVMCellCentered");

//////////////////////////////////////////////////////////////////////////////

FVMCCSparsity::FVMCCSparsity() : GlobalJacobianSparsity()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

FVMCCSparsity::~FVMCCSparsity()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCSparsity::computeNNz(std::valarray<CFint>& nnz, std::valarray<CFint>& ghostNnz)
{
  CFAUTOTRACE;

  cf_assert(socket_states.isConnected());
  cf_assert(socket_bStatesNeighbors.isConnected());

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<std::valarray<State*> > bStatesNeighbors = socket_bStatesNeighbors.getDataHandle();

  cf_assert(nnz.size() == ghostNnz.size());

  SafePtr<ConnectivityTable<CFuint> > cellFaces =
    MeshDataStack::getActive()->getConnectivity("cellFaces");

  SafePtr<MapGeoToTrsAndIdx> mapGeoToTrs =
    MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");

  SafePtr<TopologicalRegionSet> innerFaces =
    MeshDataStack::getActive()->getTrs("InnerFaces");

  const CFuint nbCells = cellFaces->nbRows();
  cf_assert(nbCells == states.size());

  // count cells on the diagonal
  nnz += 1;

  // loop on cells and add contribution from neighbor cells via the face
  for (CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    const CFuint nbFacesInCell = cellFaces->nbCols(iCell);

    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace)
    {
      // if the face is not a boundary face (=> second state is
      // not a ghost one), count the corresponding neighbor
      // NOTE that at this stage of the simulation
      // the ghost states could have been not built yet
      // (that's why you first check the nbStates of this face)
      const CFuint faceID = (*cellFaces)(iCell,iFace);
      // ignore dummy states for boundary faces
      if (!mapGeoToTrs->isBGeo(faceID))
      {
        nnz[iCell]++;

        // local index in the InnerFaces TRS
        const CFuint faceIdx = mapGeoToTrs->getIdxInTrs(faceID);
        // neighbor stateID of this face
        const CFuint neighStateID = (innerFaces->getStateID(faceIdx,1) != iCell) ?
	 innerFaces->getStateID(faceIdx,1) : innerFaces->getStateID(faceIdx,0);
        cf_assert(neighStateID < states.size());
        cf_assert(neighStateID != iCell);
        if (!states[neighStateID]->isParUpdatable()) {
          ghostNnz[iCell]++;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCSparsity::computeMatrixPattern(
  std::valarray<CFint>& nnz,
  std::valarray<CFint>& ghostNnz,
  vector< vector<CFuint> >& matrixPattern)
{
  CFAUTOTRACE;

  throw Common::NotImplementedException (FromHere(),"GlobalJacobianSparsity::computeMatrixPattern() was not implemented yet.");
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCSparsity::computeMatrixPattern
(DataSocketSink<Framework::State*, Framework::GLOBAL> statesSocket,
 Common::ConnectivityTable<CFuint>& matrixPattern)
{  
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = statesSocket.getDataHandle();
  
  SafePtr<ConnectivityTable<CFuint> > cellFaces =
    MeshDataStack::getActive()->getConnectivity("cellFaces");

  SafePtr<MapGeoToTrsAndIdx> mapGeoToTrs =
    MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");

  SafePtr<TopologicalRegionSet> innerFaces =
    MeshDataStack::getActive()->getTrs("InnerFaces");

  const CFuint nbCells = cellFaces->nbRows();
  cf_assert(nbCells == states.size());
  valarray<CFuint> nbRows(nbCells);
  nbRows = 0;
  
  // loop on cells and add contribution from neighbor cells via the face
  for (CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    const CFuint nbFacesInCell = cellFaces->nbCols(iCell);
    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace)
    {
      // if the face is not a boundary face (=> second state is
      // not a ghost one), count the corresponding neighbor
      // NOTE that at this stage of the simulation
      // the ghost states could have been not built yet
      // (that's why you first check the nbStates of this face)
      const CFuint faceID = (*cellFaces)(iCell,iFace);
      // ignore dummy states for boundary faces
      if (!mapGeoToTrs->isBGeo(faceID))
      {
	nbRows[iCell]++;
      }
    }
  }
  
  matrixPattern.resize(nbRows);
  
  // loop on cells and add contribution from neighbor cells via the face
  for (CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    const CFuint nbFacesInCell = cellFaces->nbCols(iCell);
    CFuint countF = 0;
    
    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace)
    {
      // if the face is not a boundary face (=> second state is
      // not a ghost one), count the corresponding neighbor
      // NOTE that at this stage of the simulation
      // the ghost states could have been not built yet
      // (that's why you first check the nbStates of this face)
      const CFuint faceID = (*cellFaces)(iCell,iFace);
      // ignore dummy states for boundary faces
      if (!mapGeoToTrs->isBGeo(faceID))
      {
	// local index in the InnerFaces TRS
        const CFuint faceIdx = mapGeoToTrs->getIdxInTrs(faceID);
        // neighbor stateID of this face
        const CFuint neighStateID = (innerFaces->getStateID(faceIdx,1) != iCell) ?
         innerFaces->getStateID(faceIdx,1) : innerFaces->getStateID(faceIdx,0);
        cf_assert(neighStateID < states.size());
        cf_assert(neighStateID != iCell);
	matrixPattern(iCell,countF) = neighStateID;
        countF++;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCSparsity::computeNNzNodeBased(std::valarray<CFint>& nnz, 
					std::valarray<CFint>& ghostNnz)
{
  CFAUTOTRACE;

  cf_assert(socket_states.isConnected());
  cf_assert(socket_nodes.isConnected());
  cf_assert(nnz.size() == ghostNnz.size());
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    socket_nodes.getDataHandle();
  
  std::vector<std::string> tags;
  tags.push_back ( "inner" );
  tags.push_back ( "cell" );
  
  std::vector< Common::SafePtr<TopologicalRegionSet> > innerCellsList = MeshDataStack::getActive()->getFilteredTrsList(tags);
  
  const CFuint nbNodes = nodes.size();
  vector< vector< CFuint > > neighborCells(nbNodes);
  vector< vector< Common::SafePtr<TopologicalRegionSet> > > neighborCellsTrs(nbNodes);
  
  // detect all the neighbor Cells for each Node
  const CFuint nbGroups = innerCellsList.size();
  for(CFuint iGroup=0; iGroup < nbGroups; ++iGroup){
    Common::SafePtr<TopologicalRegionSet> currTrs = innerCellsList[iGroup];
    SafePtr< ConnectivityTable< CFuint > > cellNodes =
      MeshDataStack::getActive()->getConnectivity("cellNodes_" + currTrs->getName());
    
    const CFuint nbCells = cellNodes->nbRows();
    for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
      const CFuint nbNodesInCell = cellNodes->nbCols(iCell);
      for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
        const CFuint stateID = (*cellNodes)(iCell,iNode);
        cf_assert(stateID < nbNodes);
        neighborCells[stateID].push_back(iCell);
        neighborCellsTrs[stateID].push_back(currTrs);
      }
    }
  }

  typedef set< Node*, less< Node* > > SetOfNeighbors;

  for (CFuint iNode = 0; iNode < nbNodes; ++iNode)
  {
    const CFuint nbNeighborCells = neighborCells[iNode].size();
    SetOfNeighbors* neighborNodes = new SetOfNeighbors();
    SetOfNeighbors* ghostNeighNodes = new SetOfNeighbors();

    for (CFuint iCell = 0; iCell < nbNeighborCells; ++iCell) {
      const CFuint cellTrsID = neighborCells[iNode][iCell];
      Common::SafePtr<TopologicalRegionSet> cellTrs = neighborCellsTrs[iNode][iCell];
      SafePtr< ConnectivityTable< CFuint > > cellNodes =
        MeshDataStack::getActive()->getConnectivity("cellNodes_" + cellTrs->getName());

      const CFuint nbNodesInCell = cellNodes->nbCols(cellTrsID);
      for (CFuint jNode = 0; jNode < nbNodesInCell; ++jNode) {
        // store all the neighbor nodes + the considered node itself
        Node *const currNode = nodes[(*cellNodes)(cellTrsID,jNode)];
        neighborNodes->insert(SetOfNeighbors::value_type(currNode));
        if (!currNode->isParUpdatable()) {
          // if the state is not updatable store it in ghostNeighNodes
          ghostNeighNodes->insert(SetOfNeighbors::value_type(currNode));
        }
      }
    }
    
    // number of neighbor nodes == number of non zero entries in the matrix
    const CFuint nbNeighNodes = neighborNodes->size();
    nnz[iNode] = nbNeighNodes;
    ghostNnz[iNode] = ghostNeighNodes->size();
    deletePtr( neighborNodes );
    deletePtr( ghostNeighNodes );
  }
}
  
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
