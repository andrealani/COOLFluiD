// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "samg.h"
#include "SAMGLSS/StdSetup.hh"
#include "SAMGLSS/SAMGLSSModule.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/State.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace SAMGLSS {

MethodCommandProvider< StdSetup,SAMGLSSData,SAMGLSSModule >
  stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbStates = states.size();

  // local number of states (updatable + ghost), global number of states (sum
  // of all updatable states in all processors) and number of PhysicalModel
  // equations
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbLocal = states.getLocalSize();
  const CFuint nbGlobal = states.getGlobalSize();

  // index mapping variables, where L2S is local>SAMG 0-based mapping
  std::valarray< CFuint > L2S(nbStates);
  std::valarray< bool > ghosts(nbStates);


  if (getMethodData().isParallel()) {
    // specific to parallel setup
    getMethodData().setInterface("SAMGp");

    // local iL (how COOLFluiD adresses) to samg iu/ig mapping
    CFuint iu = 0;        // next updatable samg index
    CFuint ig = nbLocal;  // ... ghost ...
    for (CFuint iL=0;iL<nbStates;++iL)
      if (states[iL]->isParUpdatable()) {
        L2S[iL] = iu++;
        ghosts[iL] = false;
      } else {
        L2S[iL] = ig++;
        ghosts[iL] = true;
      }

    // setting solver parallel parameters: get Fortran communicator (if (int)
    // MPI_Fint casting is possible) and set ghost send/receive lists
    m_sparallel_struct& y = getMethodData().getParamsParallel();

    const std::string nsp = getMethodData().getNamespace();
    MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
    
    y.icomm = static_cast< int >(MPI_Comm_c2f(comm));
    
    setParallelVector( states.getGhostSendList(),L2S,nbEqs,
      y.nshalo,y.npsnd,y.iranksnd,y.ipts,y.isndlist );
    setParallelVector( states.getGhostReceiveList(),L2S,nbEqs,
      y.nrhalo,y.nprec,y.irankrec,y.iptr,y.ireclist );

  }
  else {
    // specific to serial setup
    getMethodData().setInterface("SAMG");

    // construct dummy index mappings
    for (CFuint iL=0;iL<nbStates;++iL) {
      L2S[iL] = iL;
      ghosts[iL] = false;
    }
  }

  // setup shared data
  getMethodData().setNbStates(nbStates,nbLocal,nbGlobal);
  getMethodData().getLocalToGlobalMapping().createMapping(L2S,ghosts);

  const std::string nsp = getMethodData().getNamespace();
  MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
  
  // setup vectors
  SAMGLSSVector& sol = getMethodData().getSolVector();
  sol.create(comm,nbStates*nbEqs,nbGlobal*nbEqs,"sol");
  sol.initialize(comm,1.);

  SAMGLSSVector& rhs = getMethodData().getRhsVector();
  rhs.create(comm,nbStates*nbEqs,nbGlobal*nbEqs,"rhs");
  rhs.initialize(comm,0.);
  
  /* setup matrix: get matrix structure and change its mapping from local to
   samg (size of nnz can change, *= nbEqs) */
  std::vector< CFint > nnz(nbStates);
  std::vector< std::vector< CFuint > > nz(nbStates);
  getStructure(nnz,nz);
  expandStructureMapping(nz,nnz,nbEqs,L2S);

  SAMGLSSMatrix& mat = getMethodData().getMatrix();
  mat.createSeqBAIJ(nbEqs,nbLocal*nbEqs,nbStates*nbEqs,0,&nnz.front(),"mat");
  if (getMethodData().isPreconstruct())
    mat.createJAStructure(nz);


  /* setting solver primary parameters: this sets critical data to both serial
     and parallel runs, therefore the if(true) */
  if (true) {
    m_sprimary_struct& p = getMethodData().getParamsPrimary();
    m_sparallel_struct& y = getMethodData().getParamsParallel();
    p.nnu  = static_cast< int >(nbLocal*nbEqs);
    p.nsys = static_cast< int >(nbEqs);

    // solution approach (first digit of nsolve)
    const int napproach = p.nsolve-static_cast< int >(p.nsolve/10)*10;

    // variable-unknown and variable-point pointers (coupled systems, 4.2)
    p.ndiu = (p.nsys>1 && napproach>=2? p.nnu+y.nrhalo:1);
    p.iu = new int[p.ndiu];
    for (int i=0;i<p.ndiu;++i)
      p.iu[i] = i%p.nsys+1;

    p.ndip = (p.nsys>1 && napproach>=3? p.nnu+y.nrhalo:1);
    p.ip = new int[p.ndip];
    for (int i=0;i<p.ndip;++i)
      p.ip[i] = i/p.nsys+1;

    // unknows requiring scaling (coupled systems, 6.3)
    p.iscale = new int[p.nsys];
    for (int i=0;i<p.nsys;++i)
      p.iscale[i] = 0;  //FIXME hardcode no variables being scaled

    // adjustments for consistency
    p.p_cmplx = ((p.nsys==1)? 0.:p.p_cmplx);
  }

  /* setting solver secondary parameters: allows fine-tuning the solver
     parameters */
  if (getMethodData().useSecondaryParameters()) {
    m_ssecondary_struct& s = getMethodData().getParamsSecondary();
    int intin;
    intin = s.LEVELX;  SAMG_SET_LEVELX(&intin);
    intin = s.NPTMN;   SAMG_SET_NPTMN(&intin);
    intin = s.NCG;     SAMG_SET_NCG(&intin);
    intin = s.NWT;     SAMG_SET_NWT(&intin);
    intin = s.NTR;     SAMG_SET_NTR(&intin);
    intin = s.NRD;     SAMG_SET_NRD(&intin);
    intin = s.NRU;     SAMG_SET_NRU(&intin);
    intin = s.NRC;     SAMG_SET_NRC(&intin);
    intin = s.NP_OPT;  SAMG_SET_NP_OPT(&intin);

    double dblin;
    dblin = s.ECG;  SAMG_SET_ECG(&dblin);
    dblin = s.EWT;  SAMG_SET_EWT(&dblin);
    dblin = s.ETR;  SAMG_SET_ETR(&dblin);
  }

  /* setting other parameters of interest: this is mainly for debugging, like
     fine-tuning node numbering (avoid!) or dumping the finest-level matrix */
  if (getMethodData().isDebug()) {
    int i;
    i = 1;   SAMG_SET_REMOVE_REDUNDANCIES(&i);
    i = 1;   SAMG_SET_DO_DUMPRESTART(&i);
    // char *c ="myfilename";
    // i = 10;  SAMG_ISET_FILNAM_DUMP(static_cast< int * >(&c[0]),&i);
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::setParallelVector(
  const std::vector< std::vector< unsigned int > >& vv,
  const std::valarray< CFuint >& MAP, const CFuint nbEqs,
  int& nhalo, int& np, int*& irank, int*& ipt, int*& ilist )
{
  nhalo = 0;
  np = 0;
  for (CFuint r=0;r<vv.size();++r) {
    if (vv[r].empty()) continue;
    nhalo += vv[r].size();
    ++np;
  }
  nhalo *= static_cast< int >(nbEqs);
  irank = new int[np];
  ipt = new int[np+1];
  ilist = new int[nhalo];
  int *pir = irank;
  int *pip = ipt;
  int *pil = ilist;
  pip[0] = 1;
  for (CFuint r=0;r<vv.size();++r) {
    if (vv[r].empty()) continue;
    *pir = r;
    ++pir;
    ++pip;
    *pip = *(pip-1) + vv[r].size()*static_cast< int >(nbEqs);
    for (CFuint n=0;n<vv[r].size();++n) {
      for (CFuint e=0;e<nbEqs;++e) {
        *pil = MAP[ vv[r][n] ]*nbEqs+e +1;
        ++pil;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::getStructure(
  std::vector< CFint >& nnz,
  std::vector< std::vector< CFuint > >& nz )
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
  for (CFuint iCell=0;iCell<nb_cells;++iCell) {
    const CFuint nb_cell_states = cell_states->nbCols(iCell);
    for (CFuint iState=0;iState<nb_cell_states;++iState) {
      const State * currState = states[(*cell_states)(iCell,iState)];

      const CFuint iL = currState->getLocalID();
      for (CFuint iNeigh=0;iNeigh<nb_cell_states;++iNeigh)
        if (iNeigh!=iState) {
          const CFuint neighLocalID =
            states[(*cell_states)(iCell,iNeigh)]->getLocalID();
          nz[iL].push_back(neighLocalID);
        }
    }
  }

  /* add states to it's neighbor list (it's a non-zero entry), sort, remove
     duplicates and count the list (nnz); if state is on the boundary store
     it's neighbors local index (bStatesNeighbors) */
  for (CFuint iL=0;iL<nbStates;++iL) {
    nz[iL].push_back(iL);
    sort(nz[iL].begin(),nz[iL].end());
    nz[iL].erase(unique(nz[iL].begin(),nz[iL].end()),nz[iL].end());
    nnz[iL] = nz[iL].size();

    // store state neighbors in local indexing
    if (is_b_state[iL]) {
      bStatesNeighbors[iL].resize(nnz[iL]);
      for (int n=0;n<nnz[iL];++n)
        bStatesNeighbors[iL][n] = states[nz[iL][n]];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::expandStructureMapping(
  std::vector< std::vector< CFuint > >& nz1, std::vector< CFint >& nnz1,
  const CFuint nbEqs, const std::valarray< CFuint >& MAP )
{
  std::vector< CFint > nnz2(nnz1.size()*nbEqs);
  std::vector< std::vector< CFuint > > nz2(nz1.size()*nbEqs);

  for (CFuint iL=0;iL<nnz1.size();++iL) {   // unmapped index
    for (CFuint iEq=0;iEq<nbEqs;++iEq) {
      const CFuint iS = MAP[iL]*nbEqs+iEq;  // mapped index

      nnz2[iS] = nnz1[iL]*nbEqs;
      nz2[iS].resize(nnz1[iL]*nbEqs);
      for (int iN=0;iN<nnz1[iL];++iN) {  // cycle neighbor list
        const CFuint neigh_mapped = MAP[nz1[iL][iN]]*nbEqs;
        for (CFuint iEq2=0;iEq2<nbEqs;++iEq2)
          nz2[iS][iN*nbEqs+iEq2] = neigh_mapped+iEq2;
      }
      sort(nz2[iS].begin(),nz2[iS].end());
    }
  }
  nnz1.swap(nnz2);
  nz1.swap(nz2);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SAMGLSS
}  // namespace COOLFluiD

