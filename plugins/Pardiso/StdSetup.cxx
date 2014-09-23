// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Pardiso/StdSetup.hh"
#include "Pardiso/PardisoModule.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/State.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace Pardiso {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSetup,PardisoData,PardisoModule >
  stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;

  PardisoData& d = getMethodData();

  // number of states
  DataHandle < Framework::State*, Framework::GLOBAL > states =
    socket_states.getDataHandle();
  const CFuint nbStates = states.size();

  // number of variables per state
  CFuint nbEqs = 0;
  Common::SafePtr< std::valarray< bool > > maskarray = d.getMaskArray();
  for (CFuint i=0; i<maskarray->size(); ++i)
    if ((*maskarray)[i])
      ++nbEqs;

  // setup local to internal index mapping (no renumbering is performed)
  //TODO: maybe use createStdSequential()?
  std::valarray< CFuint > L2I(nbStates);
  std::valarray< bool > ghosts(nbStates);
  for (CFuint iL=0; iL<nbStates; ++iL) {
    L2I[iL] = iL;
    ghosts[iL] = false;
  }
  d.setNbStates(nbStates);
  d.setNbEquations(nbEqs);
  d.getLocalToGlobalMapping().createMapping(L2I,ghosts);

  // setup system vectors
#ifndef CF_HAVE_MPI
  // to compile without MPI (as seen in from Framework/LSSVector.hh)  
  const MPI_Comm MPI_COMM_WORLD = 0;
#endif
  PardisoVector& sol = d.getSolVector();
  PardisoVector& rhs = d.getRhsVector();
  sol.create(MPI_COMM_WORLD,nbStates*nbEqs,nbStates*nbEqs,"sol");
  rhs.create(MPI_COMM_WORLD,nbStates*nbEqs,nbStates*nbEqs,"rhs");
  sol.initialize(MPI_COMM_WORLD,0.);
  rhs.initialize(MPI_COMM_WORLD,0.);

  // setup system matrix
  std::vector< CFint > nnz(nbStates);
  std::vector< std::vector< CFuint > > nz(nbStates);
  getStructure(nnz,nz);
  expandStructureMapping(nz,nnz,nbEqs);
  PardisoMatrix& mat = d.getMatrix();
  mat.createSeqBAIJ(nbEqs,nbStates*nbEqs,nbStates*nbEqs,0,&nnz.front(),"mat");
  mat.createJAStructure(nz);
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::getStructure(
  std::vector< CFint >& nnz, std::vector< std::vector< CFuint > >& nz )
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<std::valarray< State* > > bStatesNeighbors =
    socket_bStatesNeighbors.getDataHandle();

  const CFuint nbStates = states.size();
  cf_assert(nz.size()==nbStates);
  cf_assert(nnz.size()==nbStates);

  // loop over all the boundary TRSs to detect all the boundary states
  std::valarray< bool > is_b_state(false,nbStates);
  std::vector< Common::SafePtr<TopologicalRegionSet> > alltrs =
    MeshDataStack::getActive()->getTrsList();
  std::vector< Common::SafePtr<TopologicalRegionSet> >::const_iterator itrs;
  for (itrs = alltrs.begin(); itrs != alltrs.end(); ++itrs) {
    Common::SafePtr<TopologicalRegionSet> currTrs = *itrs;
    if (currTrs->getName()=="InnerCells" || currTrs->getName()=="InnerFaces")
      continue;
    Common::SafePtr< std::vector< CFuint > > bStates =
      currTrs->getStatesInTrs();
    for (CFuint i = 0; i < bStates->size(); ++i)
      is_b_state[(*bStates)[i]] = true;
  }

  // store the states from all cells
  Common::SafePtr< Common::ConnectivityTable< CFuint > > cell_states =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");
  const CFuint nb_cells = cell_states->nbRows();
  for (CFuint iCell=0; iCell<nb_cells; ++iCell) {
    const CFuint nb_cell_states = cell_states->nbCols(iCell);
    for (CFuint iState=0; iState<nb_cell_states; ++iState) {
      const State * currState = states[(*cell_states)(iCell,iState)];

      const CFuint iL = currState->getLocalID();
      for (CFuint iNeigh=0; iNeigh<nb_cell_states; ++iNeigh)
        if (iNeigh!=iState) {
          const CFuint neighLocalID =
            states[(*cell_states)(iCell,iNeigh)]->getLocalID();
          nz[iL].push_back(neighLocalID);
        }
    }
  }

  // add states to it's neighbor list (it's a non-zero entry), sort, remove
  // duplicates and count the list (nnz); if state is on the boundary store
  // it's neighbors local index (bStatesNeighbors)
  for (CFuint iL=0; iL<nbStates; ++iL) {
    nz[iL].push_back(iL);
    sort(nz[iL].begin(),nz[iL].end());
    nz[iL].erase(unique(nz[iL].begin(),nz[iL].end()),nz[iL].end());
    nnz[iL] = nz[iL].size();

    // store state neighbors in local indexing
    if (is_b_state[iL]) {
      bStatesNeighbors[iL].resize(nnz[iL]);
      for (int n=0; n<nnz[iL]; ++n)
        bStatesNeighbors[iL][n] = states[nz[iL][n]];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::expandStructureMapping(
  std::vector< std::vector< CFuint > >& nz1, std::vector< CFint >& nnz1,
  const CFuint nbEqs )
{
  std::vector< CFint > nnz2(nnz1.size()*nbEqs);
  std::vector< std::vector< CFuint > > nz2(nz1.size()*nbEqs);

  for (CFuint iL=0; iL<nnz1.size(); ++iL) {   // unmapped index
    for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
      const CFuint iS = iL*nbEqs+iEq;         // mapped index

      nnz2[iS] = nnz1[iL]*nbEqs;
      nz2[iS].resize(nnz1[iL]*nbEqs);
      for (int iN=0; iN<nnz1[iL]; ++iN) {     // cycle neighbor list
        const CFuint neigh_mapped = nz1[iL][iN]*nbEqs;
        for (CFuint iEq2=0; iEq2<nbEqs; ++iEq2)
          nz2[iS][iN*nbEqs+iEq2] = neigh_mapped+iEq2;
      }
      sort(nz2[iS].begin(),nz2[iS].end());
    }
  }
  nnz1.swap(nnz2);
  nz1.swap(nz2);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Pardiso
} // namespace COOLFluiD

