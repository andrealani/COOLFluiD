#include "Environment/ObjectProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCCSparsityNoBlock.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FVMCCSparsityNoBlock,
               GlobalJacobianSparsity,
               FiniteVolumeModule>
aFVMCCSparsityNoBlockProvider("FVMCellCenteredNoBlock");

//////////////////////////////////////////////////////////////////////////////

FVMCCSparsityNoBlock::FVMCCSparsityNoBlock() : GlobalJacobianSparsity()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

FVMCCSparsityNoBlock::~FVMCCSparsityNoBlock()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCSparsityNoBlock::computeNNz(std::valarray<CFint>& nnzOut, std::valarray<CFint>& ghostNnzOut)
{
  CFAUTOTRACE;

  cf_assert(socket_states.isConnected());
  cf_assert(socket_bStatesNeighbors.isConnected());

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<std::valarray<State*> > bStatesNeighbors = socket_bStatesNeighbors.getDataHandle();

  cf_assert(nnzOut.size() == ghostNnzOut.size());

  SafePtr<ConnectivityTable<CFuint> > cellFaces =
    MeshDataStack::getActive()->getConnectivity("cellFaces");

  SafePtr<MapGeoToTrsAndIdx> mapGeoToTrs =
    MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");

  SafePtr<TopologicalRegionSet> innerFaces =
    MeshDataStack::getActive()->getTrs("InnerFaces");

  const CFuint nbCells = cellFaces->nbRows();
  cf_assert(nbCells == states.size());

  // count cells on the diagonal
  valarray<CFint> nnz(nnzOut.size());
  valarray<CFint> ghostNnz(nnzOut.size());
  
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

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint outSize = nnz.size()*nbEqs;
  nnzOut.resize(outSize);
  ghostNnzOut.resize(outSize);
  
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    const CFuint startRow = iCell*nbEqs;
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      nnzOut[startRow+iEq] = nnz[iCell]*nbEqs;
      ghostNnzOut[startRow+iEq] = ghostNnz[iCell]*nbEqs;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCSparsityNoBlock::computeMatrixPattern(
  std::valarray<CFint>& nnz,
  std::valarray<CFint>& ghostNnz,
  vector< vector<CFuint> >& matrixPattern)
{
  CFAUTOTRACE;

  throw Common::NotImplementedException (FromHere(),"FVMCCSparsityNoBlock::computeMatrixPattern() was not implemented yet.");
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCSparsityNoBlock::computeMatrixPattern
(DataSocketSink<Framework::State*, Framework::GLOBAL> statesSocket,
 Common::ConnectivityTable<CFuint>& matrixPattern)
{  
  CFAUTOTRACE;
  throw Common::NotImplementedException (FromHere(),"GlobalJacobianSparsity::computeMatrixPattern() was not implemented yet.");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
