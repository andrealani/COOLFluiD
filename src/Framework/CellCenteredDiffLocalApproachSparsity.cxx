// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include "Common/NotImplementedException.hh"

#include "Environment/ObjectProvider.hh"

#include "Framework/CellCenteredDiffLocalApproachSparsity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/State.hh"
#include "Framework/CFSide.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

 namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CellCenteredDiffLocalApproachSparsity,
               GlobalJacobianSparsity,
               FrameworkLib>
aCellCenteredDiffLocalApproachSparsityProvider("CellCenteredDiffLocalApproach");

//////////////////////////////////////////////////////////////////////////////

CellCenteredDiffLocalApproachSparsity::CellCenteredDiffLocalApproachSparsity() : GlobalJacobianSparsity()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

CellCenteredDiffLocalApproachSparsity::~CellCenteredDiffLocalApproachSparsity()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void CellCenteredDiffLocalApproachSparsity::computeNNz(std::valarray<CFint>& nnz,
                                      std::valarray<CFint>& ghostNnz)
{
  CFAUTOTRACE;

  cf_assert(nnz.size() == ghostNnz.size());

  cf_assert(socket_states.isConnected());

  // get states datahandle
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // get InnerCells and InnerFaces TRSs
  SafePtr<TopologicalRegionSet> innerCells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  SafePtr<TopologicalRegionSet> innerFaces =
    MeshDataStack::getActive()->getTrs("InnerFaces");

  // get connectivities
  SafePtr< ConnectivityTable<CFuint> > faceToCells =
    MeshDataStack::getActive()->getConnectivity("InnerFaces-Faces2Cells");

  SafePtr< ConnectivityTable< CFuint > > cellToFaces =
      MeshDataStack::getActive()->getConnectivity("InnerCells-Cells2Faces");

  SafePtr< ConnectivityTable< CFuint > > cellStates =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  SafePtr< ConnectivityTable< CFuint > > faceIntFaceIdx =
      MeshDataStack::getActive()->getConnectivity("faceToInFaceIdxOrientOrBCIdx");

  // number of cells
  const CFuint nbCells = cellStates->nbRows();
  cf_assert(nbCells == innerCells->getLocalNbGeoEnts());

  // variable to keep track of cells that are a part of the numerical scheme for a cell
  vector< vector< CFuint > > cellsInCellScheme(nbCells);

  // loop on cells and add contributions
  for (CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    // number of states in this cell
    const CFuint nbStates = cellStates->nbCols(iCell);
    cf_assert(nbStates > 0);

    // add cell to numerical scheme
    cellsInCellScheme[iCell].push_back(iCell);

    // add number of non-zero contributions to each state in cell
    const CFuint firstStateID = (*cellStates)(iCell,0);
    if (states[firstStateID]->isParUpdatable()) // if the first cell-state is parUpdatable, they all should be
    {
      for (CFuint iState = 0; iState < nbStates; ++iState)
      {
        const CFuint stateID = (*cellStates)(iCell,iState);
        cf_assert(stateID < states.size());
        cf_assert(states[stateID]->isParUpdatable());
        cf_assert(stateID < nnz.size());
        nnz[stateID] += nbStates;
      }
    }
    else
    {
      for (CFuint iState = 0; iState < nbStates; ++iState)
      {
        const CFuint stateID = (*cellStates)(iCell,iState);
        cf_assert(stateID < states.size());
        cf_assert(!states[stateID]->isParUpdatable());
        cf_assert(stateID < ghostNnz.size());
        ghostNnz[stateID] += nbStates;
      }
    }
  }

  // number of inner faces
  const CFuint nbInnerFaces = faceToCells->nbRows();
  cf_assert(innerFaces->getLocalNbGeoEnts() == nbInnerFaces);

  // loop on faces and add contributions of face neighbours
  for (CFuint iFace = 0; iFace < nbInnerFaces; ++iFace)
  {
    // neighbouring cell IDs and number of states in cells
    vector< CFuint > cellIDs(2);
    vector< CFuint > nbStatesCells(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      cellIDs[iSide] = (*faceToCells)(iFace,iSide);
      cf_assert(cellIDs[iSide] < nbCells);
      nbStatesCells[iSide] = cellStates->nbCols(cellIDs[iSide]);
    }

    // add neigbouring cell contributions
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      // variable for other side
      const CFuint iOtherSide = iSide == LEFT ? RIGHT : LEFT;

      // add cells to numerical scheme
      cellsInCellScheme[cellIDs[iSide]].push_back(cellIDs[iOtherSide]);

      // add contribution of states in neighbouring cell
      const CFuint firstStateID = (*cellStates)(cellIDs[iOtherSide],0);
      if (states[firstStateID]->isParUpdatable())
      // if the first cell-state is parUpdatable, then all should be!!
      {
        for (CFuint iState = 0 ; iState < nbStatesCells[iSide]; ++iState)
        {
          const CFuint stateID = (*cellStates)(cellIDs[iSide],iState);
          cf_assert(stateID < states.size());
          cf_assert(stateID < nnz.size());
          nnz[stateID] += nbStatesCells[iOtherSide];
        }
      }
      else
      // if the first cell-state is ghost, then all should be!!
      {
        for (CFuint iState = 0 ; iState < nbStatesCells[iSide]; ++iState)
        {
          const CFuint stateID = (*cellStates)(cellIDs[iSide],iState);
          cf_assert(stateID < states.size());
          cf_assert(stateID < ghostNnz.size());
          ghostNnz[stateID] += nbStatesCells[iOtherSide];
        }
      }
    }
  } // loop inner faces

  // loop on faces and add contributions of the neighbours of the face neighbours
  for (CFuint iFace = 0; iFace < nbInnerFaces; ++iFace)
  {
    // neighbouring cell IDs and number of states in cells
    vector< CFuint > cellIDs(2);
    vector< CFuint > nbStatesCells(2);
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      cellIDs[iSide] = (*faceToCells)(iFace,iSide);
      cf_assert(cellIDs[iSide] < nbCells);
      nbStatesCells[iSide] = cellStates->nbCols(cellIDs[iSide]);
    }

    // add neigbouring cell contributions
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      // variable for other side
      const CFuint iOtherSide = iSide == LEFT ? RIGHT : LEFT;

      // loop over neighbouring faces
      const CFuint nbrNeighbFaces = cellToFaces->nbCols(cellIDs[iSide]);
      for (CFuint iNFace = 0; iNFace < nbrNeighbFaces; ++iNFace)
      {
        // get face ID
        const CFuint faceID = (*cellToFaces)(cellIDs[iSide],iNFace);

        // check if face is inner face
        const bool isInnerFace = faceIntFaceIdx->nbCols(faceID) == 2;

        // if inner face, check the neighbouring cell
        if (isInnerFace)
        {
          // get inner face index
          const CFuint intFaceIdx = (*faceIntFaceIdx)(faceID,0);

          // check if face is not current inner face
          const bool isOtherNeighbourFace = intFaceIdx != iFace;

          // if other neighbouring face, add neighbouring cell to numerical scheme
          if (isOtherNeighbourFace)
          {
            // get inner face neighbouring cell other than current cell
            cf_assert(intFaceIdx < faceToCells->nbRows());
            const CFuint neighbCellID =
                (*faceToCells)(intFaceIdx,LEFT ) == cellIDs[iSide] ?
                (*faceToCells)(intFaceIdx,RIGHT) :
                (*faceToCells)(intFaceIdx,LEFT );
            cf_assert(neighbCellID < cellStates->nbRows());

            // check if cell is not in numerical scheme yet
            bool isNotInNumScheme = true;
            const CFuint nbrCellsInScheme = cellsInCellScheme[cellIDs[iOtherSide]].size();
            for (CFuint iCell = 0; iCell < nbrCellsInScheme && isNotInNumScheme; ++iCell)
            {
              isNotInNumScheme =
                  cellsInCellScheme[cellIDs[iOtherSide]][iCell] != neighbCellID;
            }

            // if cell is not yet in numerical scheme, add the cell
            if (isNotInNumScheme)
            {
              // add to list
              cellsInCellScheme[cellIDs[iOtherSide]].push_back(neighbCellID);

              // get number of states in the cell
              const CFuint nbrNeighbCellStates = cellStates->nbCols(neighbCellID);

              // add nonzero entries to states
              const CFuint firstStateID = (*cellStates)(neighbCellID,0);
              if (states[firstStateID]->isParUpdatable())
              // if the first cell-state is parUpdatable, then all should be!!
              {
                for (CFuint iState = 0 ; iState < nbStatesCells[iOtherSide]; ++iState)
                {
                  const CFuint stateID = (*cellStates)(cellIDs[iOtherSide],iState);
                  cf_assert(stateID < states.size());
                  cf_assert(stateID < nnz.size());
                  nnz[stateID] += nbrNeighbCellStates;
                }
              }
              else
              // if the first cell-state is ghost, then all should be!!
              {
                for (CFuint iState = 0 ; iState < nbStatesCells[iOtherSide]; ++iState)
                {
                  const CFuint stateID = (*cellStates)(cellIDs[iOtherSide],iState);
                  cf_assert(stateID < states.size());
                  cf_assert(stateID < ghostNnz.size());
                  ghostNnz[stateID] += nbrNeighbCellStates;
                }
              }
            }
          }
        }
      }
    }
  } // loop inner faces

//   for (CFuint iCell = 0; iCell < nbCells; ++iCell)
//   {
//     CF_DEBUG_OBJ(iCell);
//     CF_DEBUG_OBJ(nnz[iCell]);
//     CF_DEBUG_OBJ(ghostNnz[iCell]);
// /*    for (CFuint i = 0; i < cellsInCellScheme[iCell].size(); ++i)
//     {
//       CF_DEBUG_OBJ(cellsInCellScheme[iCell][i]);
//     }*/
//   }

  ofstream fout ("nnz.dat");

        fout << " NNZ " << "\n";
      for (CFuint i = 0 ; i < nnz.size(); ++i)
      {
        fout << nnz[i] << "\n";
      }
        fout << "Ghost NNZ " << "\n";
      for (CFuint i = 0 ; i < ghostNnz.size(); ++i)
      {
        fout << ghostNnz[i] << "\n";
      }


}

//////////////////////////////////////////////////////////////////////////////

void CellCenteredDiffLocalApproachSparsity::computeMatrixPattern(
  std::valarray<CFint>& nnz,
  std::valarray<CFint>& ghostNnz,
  vector< vector<CFuint> >& matrixPattern)
{
  CFAUTOTRACE;
  throw Common::NotImplementedException (FromHere(),"GlobalJacobianSparsity::computeMatrixPattern() was not implemented yet.");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
