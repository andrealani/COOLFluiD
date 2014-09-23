#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "PLaS/PLaSModule.hh"
#include "PLaS/PLaSTrackingData.hh"
#include "PLaS/StgImplementationBlockable.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider< StgImplementationBlockable,PLaSTrackingData,StgImplementation,PLaSModule > plasImplementationBlockableProvider("StgImplementationBlockable");

//////////////////////////////////////////////////////////////////////////////

StgImplementationBlockable::StgImplementationBlockable(const std::string& name) :
    StgImplementation(name),
    s_nodes("nodes"),                  // socket sinks
    s_states("states"),                // ...
    s_faceneighcell("faceNeighCell"),  // ...
    s_sol_tm0("SolutionAt(t-0)"),      // ...
    s_ddx_tm0("SolutionAt(t-0)ddx"),   // ...
    s_ddy_tm0("SolutionAt(t-0)ddy"),   // ...
    s_ddz_tm0("SolutionAt(t-0)ddz"),   // ...
    s_sol_tm1("SolutionAt(t-1)"),      // ...
    s_ddx_tm1("SolutionAt(t-1)ddx"),   // ...
    s_ddy_tm1("SolutionAt(t-1)ddy"),   // ...
    s_ddz_tm1("SolutionAt(t-1)ddz"),   // ...
    s_nvolume("NodalVolume"),          // ...
    s_evolume("ElementVolumes"),       // ...
    s_ienormals("ElementNormals"),     // ...
    s_benormals("BoundariesNormals")   // ...
{
}

//////////////////////////////////////////////////////////////////////////////

void StgImplementationBlockable::setup()
{
  StgImplementation::setup();
  getMethodData().m_blockable = true;
}

//////////////////////////////////////////////////////////////////////////////

void StgImplementationBlockable::setupVolumes(int& nelm, double& vtot, double& vmin, double& vmax)
{
  CFAUTOTRACE;

  DataHandle< CFreal > h_elementvolumes = s_evolume.getDataHandle();
  nelm = (int) h_elementvolumes.size();
  vtot = 0.;
  vmin = 1.e99;
  vmax = 0.;
  for (int ic=0; ic<nelm; ++ic) {
    const double v = (double) h_elementvolumes[ic];
    vtot += v;
    vmin = std::min(vmin,v);
    vmax = std::max(vmax,v);
  }
}

//////////////////////////////////////////////////////////////////////////////

SafePtr< ConnectivityTable< CFuint > > StgImplementationBlockable::setupConnectivity()
{
  // return connectivity of first trs found
  std::vector< SafePtr< TopologicalRegionSet > > inner =
    MeshDataStack::getActive()->getFilteredTrsList("inner");
  cf_assert_desc("no trs with inner tag found!",inner.size());
  return inner[0]->getGeo2NodesConn();
}

//////////////////////////////////////////////////////////////////////////////

void StgImplementationBlockable::setupElementToElement()
{
  CFAUTOTRACE;
  // this sets element to element connectivity, for those sharing faces. can
  // be extended easily to sharing edges or nodes, but PLaS wants it this way.
  // for non-existing neighbors sets -1.

  // set number of element neighbors (boundaries require special treatment)
  const CFuint Nelm = m_geo2nodes->nbRows();
  const CFuint Nnod = s_nodes.getDataHandle().size();
  const CFuint nbFaces = (CFuint) nbDim+1;       // simplexes have from 1 to nbDim+1 neighbors
  std::valarray< CFuint > nneigh(nbFaces,Nelm);  // size of connectivity table

  // set nodes to elements connectivity
  std::vector< std::vector< CFuint > >
    nodes2geo(Nnod,std::vector< CFuint >());
  for (CFuint ic=0; ic<Nelm; ++ic)
    for (CFuint in=0; in<nbFaces; ++in)
      nodes2geo[ (*m_geo2nodes)(ic,in) ].push_back(ic);
  for (CFuint in=0; in<Nnod; ++in) {
    sort(nodes2geo[in].begin(),nodes2geo[in].end());
    nodes2geo[in].erase(
      unique( nodes2geo[in].begin(),nodes2geo[in].end() ),
      nodes2geo[in].end() );
  }

  // set neighbors entries
  m_geo2geo.resize(nneigh,-1);
  for (CFuint ic=0; ic<Nelm; ++ic) {

    // get all touching elements (copy elements indices and remove duplicates)
    std::vector< CFuint > ngeos;
    for (CFuint in=0; in<nbFaces; ++in)
      ngeos.insert(
        ngeos.end(),
        nodes2geo[(*m_geo2nodes)(ic,in)].begin(),
        nodes2geo[(*m_geo2nodes)(ic,in)].end() );
    sort(ngeos.begin(),ngeos.end());
    ngeos.erase( unique(ngeos.begin(),ngeos.end()), ngeos.end() );

    // keep track of shared nodes of touching elements
    std::vector< std::vector< bool > >
      nshared(ngeos.size(),std::vector< bool >(nbFaces,false));
    for (CFuint in=0; in<nbFaces; ++in)
      for (CFuint ic2=0; ic2<ngeos.size(); ++ic2) {
        std::vector< CFuint > nn(nbFaces,0);
        for (CFuint in2=0; in2<nbFaces; ++in2)
          nn[in2] = (*m_geo2nodes)(ngeos[ic2],in2);
        nshared[ic2][in] = count(nn.begin(),nn.end(),(*m_geo2nodes)(ic,in))>0;
      }

    // for each touching element, get ones sharing exactly nbDim nodes (they
    // should share a face). index f is where it is stored in the table, which
    // is the face index shared by the elements according to gambit convention
    // (which matches coolfluid, by the way.)
    for (CFuint ic2=0; ic2<ngeos.size(); ++ic2)
      if (count(nshared[ic2].begin(),nshared[ic2].end(),true)==nbDim) {
        // the next comments are the gambit convention; for simplices and in
        // this case PLaS uses "face opposing node", otherwise follows gambit.
        const int f = (nbDim>2?
          (!nshared[ic2][0]? /*2*/ 0 :
          (!nshared[ic2][1]? /*3*/ 1 :
          (!nshared[ic2][2]? /*1*/ 2 :
          (!nshared[ic2][3]? /*0*/ 3 : -1 )))) :
          (!nshared[ic2][0]? /*1*/ 0 :
          (!nshared[ic2][1]? /*2*/ 1 :
          (!nshared[ic2][2]? /*0*/ 2 : -1 ))) );
        cf_assert_desc("element neighbors: impossible face!",f>=0);
        cf_assert_desc("element neighbors: overwriting!",m_geo2geo(ic,f)<0);
        m_geo2geo(ic,f) = ngeos[ic2];
        --nneigh[ic];
      }

  }
}

//////////////////////////////////////////////////////////////////////////////

void StgImplementationBlockable::PLaS_SetFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp)
{
  CFAUTOTRACE;

  // set number of dimensions and variables
  nbDim = (int) s_nodes.getDataHandle()[0]->size();
  nbVar = s_states.getDataHandle()[0]->size();

  // setup element to node connectivity table
  m_geo2nodes = setupConnectivity();
  cf_assert_desc(
    "problem setting connectivity table (tag:\"inner\")!",
    m_geo2nodes.isNotNull() );

  // setup element to element connectivity table
  setupElementToElement();

  // setup number of elements and volume statistics
  setupVolumes(
    fp->numElm,          // (int) number of elements
    fp->domainVolume,    // (double) total volume of domain
    fp->minElmVolume,    // (double) smallest element volume of domain
    fp->maxElmVolume );  // (double) largest element volume of domain
}

//////////////////////////////////////////////////////////////////////////////

void StgImplementationBlockable::PLaS_SetFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->time = (double) SubSystemStatusStack::getActive()->getCurrentTime();  // (double) current time
  fp->iter = (int)    SubSystemStatusStack::getActive()->getNbIter();       // (int) current iteration
  fp->writeOutput = 1;  // (int) flag whether the flow solver writes outputs
}

//////////////////////////////////////////////////////////////////////////////

void StgImplementationBlockable::PLaS_SetPartitioningData(PLAS_PART_DATA *part)
{
  PLaSTrackingData& d = getMethodData();

  // (int) number of updatable/ghost node pairs
  part->numNodPairs = d.m_ghosts_srank.size();
  if (!part->numNodPairs) {
    part->sendRank    = CFNULL;  // (int*) sending ranks
    part->sendNodeIdx = CFNULL;  // (int*) ... rank local node numbers
    part->recvRank    = CFNULL;  // (int*) receiving ranks
    part->recvNodeIdx = CFNULL;  // (int*) ... rank local node numbers
  }
  else {
    part->sendRank    = &d.m_ghosts_srank[0];
    part->sendNodeIdx = &d.m_ghosts_sindex[0];
    part->recvRank    = &d.m_ghosts_rrank[0];
    part->recvNodeIdx = &d.m_ghosts_rindex[0];
  }
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationBlockable::PLaS_GetElmNode(int elm, int enod)
{
  return (*m_geo2nodes)((CFuint) elm,(CFuint) enod);
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationBlockable::PLaS_GetElmNeighbour(int elm, int eface)
{
  return m_geo2geo(elm,eface);
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetBndFaceRefCoord(int bnd, int bface, int dim)
{
  return PLaS_GetNodCoord(
    (int) (*getMethodData().m_boundary[bnd]->getGeo2NodesConn())(bface,0),
    dim );
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationBlockable::PLaS_GetBndDomElm(int bnd, int bface, int dummy)
{
  return (int) s_faceneighcell.getDataHandle()[
    getMethodData().m_boundary[bnd]->getLocalGeoID((CFuint) bface) ].first;
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetNodCoord(int nod, int dim)
{
  return (*(s_nodes.getDataHandle())[nod])[dim];
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetElmNormComp(int elm, int eface, int dim)
{
  return s_ienormals.getDataHandle()[elm][dim][eface];
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetBndFaceNormComp(int bnd, int bface, int dim)
{
  return s_benormals.getDataHandle()[bnd][dim][bface];
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetElmFaceMiddlePoint(int elm, int eface, int dim)
{
  // attention: this does not follow gambit convention! it is supposed to be
  // asking the face 'eface' which is the one not (!) containing element node
  // index 'eface'!
  const int Nelmnode = nbDim+1;
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  if (nbDim>2)
    return (1./3.)*(
      (*h_nodes[(*m_geo2nodes)(elm, (eface+1)%Nelmnode )])[dim] +
      (*h_nodes[(*m_geo2nodes)(elm, (eface+2)%Nelmnode )])[dim] +
      (*h_nodes[(*m_geo2nodes)(elm, (eface+3)%Nelmnode )])[dim] );
  return 0.5*(
    (*h_nodes[(*m_geo2nodes)(elm, (eface+1)%Nelmnode )])[dim] +
    (*h_nodes[(*m_geo2nodes)(elm, (eface+2)%Nelmnode )])[dim] );
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetNodVol(int nod)
{
  return s_nvolume.getDataHandle()[nod];
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetElmVol(int elm)
{
  return s_evolume.getDataHandle()[elm];
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationBlockable::PLaS_GetNumBndFaces(int bnd)
{
  return (int) getMethodData().m_boundary[bnd]->getGeo2NodesConn()->nbRows();
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationBlockable::PLaS_GetWallBndFlag(int bnd)
{
  return (getMethodData().m_boundary_is_wall[bnd]? 1:0);
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationBlockable::PLaS_GetPerBndFlag(int bnd)
{
  return 0;  // no periodic boundaries supported, no intention to ever do it
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetPerBndOffset(int bnd, int dim)
{
  return 0.;  // no periodic boundaries supported, no intention to ever do it
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationBlockable::PLaS_GetElementType(int elm)
{
  return ELM_SIMPLEX;
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetVelocityComp(int nod, int dim)
{
  if (dim<0 || dim>=nbDim)
    return 0.;
  const int i = getMethodData().m_vidx[dim];
  return (i<0? getMethodData().m_vdef[dim] : s_sol_tm0.getDataHandle()(nod,i,nbVar) );
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetVelocityCompOld(int nod, int dim)
{
  if (dim<0 || dim>=nbDim)
    return 0.;
  const int i = getMethodData().m_vidx[dim];
  return (i<0? getMethodData().m_vdef[dim] : s_sol_tm1.getDataHandle()(nod,i,nbVar) );
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetVelocityDerivativeComp(int nod, int idim, int jdim)
{
  if (idim<0 || idim>=nbDim || jdim<0 || jdim>=nbDim)
    return 0.;
  const int i = getMethodData().m_vidx[jdim];
  return (i<0? 0. :
         (idim==0? s_ddx_tm0.getDataHandle()(nod,i,nbVar) :
         (idim==1? s_ddy_tm0.getDataHandle()(nod,i,nbVar) :
                   s_ddz_tm0.getDataHandle()(nod,i,nbVar) )));
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetVelocityDerivativeCompOld(int nod, int idim, int jdim)
{
  if (idim<0 || idim>=nbDim || jdim<0 || jdim>=nbDim)
    return 0.;
  const int i = getMethodData().m_vidx[jdim];
  return (i<0? 0. :
         (idim==0? s_ddx_tm1.getDataHandle()(nod,i,nbVar) :
         (idim==1? s_ddy_tm1.getDataHandle()(nod,i,nbVar) :
                   s_ddz_tm1.getDataHandle()(nod,i,nbVar) )));
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetTemperature(int nod)
{
  const int i = getMethodData().m_tidx;
  return (i<0? getMethodData().m_tdef : s_sol_tm0.getDataHandle()(nod,i,nbVar));
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetTemperatureOld(int nod)
{
  const int i = getMethodData().m_tidx;
  return (i<0? getMethodData().m_tdef : s_sol_tm1.getDataHandle()(nod,i,nbVar));
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetPressure(int nod)
{
  const int i = getMethodData().m_pidx;
  return (i<0? getMethodData().m_pdef : s_sol_tm0.getDataHandle()(nod,i,nbVar));
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetPressureOld(int nod)
{
  const int i = getMethodData().m_pidx;
  return (i<0? getMethodData().m_pdef : s_sol_tm1.getDataHandle()(nod,i,nbVar));
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationBlockable::PLaS_StartElementSearch(double *pos)
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationBlockable::PLaS_EndElementSearch(double *pos)
{
  return (int) m_geo2nodes->nbRows()-1;
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationBlockable::PLaS_GetEulerianTimeScale(int nod)
{
  DataHandle< CFreal > h_sol = s_sol_tm0.getDataHandle();
  DataHandle< CFreal > h_ddx = s_ddx_tm0.getDataHandle();
  DataHandle< CFreal > h_ddy = s_ddy_tm0.getDataHandle();
  DataHandle< CFreal > h_ddz = s_ddz_tm0.getDataHandle();
  const int iu = getMethodData().m_vidx[0];
  const int iv = getMethodData().m_vidx[1];
  const int iw = (nbDim>2? getMethodData().m_vidx[2]:0);
  if (iu<0 || iv<0 || iw<0)
    return 1.;

  const CFreal nu = getMethodData().m_plasdata.fp.nuCont;
  const CFreal v2 =
      h_sol(nod,iu,nbVar) * h_sol(nod,iu,nbVar)
    + h_sol(nod,iv,nbVar) * h_sol(nod,iv,nbVar)
    + (nbDim>2? h_sol(nod,iw,nbVar) * h_sol(nod,iw,nbVar) : 0.);
  const CFreal ke = 0.5*v2;

  CFreal dissip = 0.;
#define UU nod,iu,nbVar
#define VV nod,iv,nbVar
#define WW nod,iw,nbVar
  if (nbDim>2) {
    dissip += (h_ddx(UU) + h_ddx(UU))*(h_ddx(UU) + h_ddx(UU)) +
              (h_ddx(VV) + h_ddy(UU))*(h_ddx(VV) + h_ddy(UU)) +
              (h_ddx(WW) + h_ddz(UU))*(h_ddx(WW) + h_ddz(UU)) +

              (h_ddy(UU) + h_ddx(VV))*(h_ddy(UU) + h_ddx(VV)) +
              (h_ddy(VV) + h_ddy(VV))*(h_ddy(VV) + h_ddy(VV)) +
              (h_ddy(WW) + h_ddz(VV))*(h_ddy(WW) + h_ddz(VV)) +

              (h_ddz(UU) + h_ddx(WW))*(h_ddz(UU) + h_ddx(WW)) +
              (h_ddz(VV) + h_ddy(WW))*(h_ddz(VV) + h_ddy(WW)) +
              (h_ddz(WW) + h_ddz(WW))*(h_ddz(WW) + h_ddz(WW));
  }
  else {
    dissip += (h_ddx(UU) + h_ddx(UU))*(h_ddx(UU) + h_ddx(UU)) +
              (h_ddx(VV) + h_ddy(UU))*(h_ddx(VV) + h_ddy(UU)) +

              (h_ddy(UU) + h_ddx(VV))*(h_ddy(UU) + h_ddx(VV)) +
              (h_ddy(VV) + h_ddy(VV))*(h_ddy(VV) + h_ddy(VV));
  }
#undef UU
#undef VV
#undef WW
  dissip *= 0.5*nu;

  return (ke/dissip);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

