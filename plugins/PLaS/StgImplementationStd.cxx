#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "PLaS/PLaSModule.hh"
#include "PLaS/PLaSTrackingData.hh"
#include "PLaS/StgImplementationStd.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider< StgImplementationStd,PLaSTrackingData,StgImplementation,PLaSModule > plasImplementationStdProvider("StgImplementationStd");

//////////////////////////////////////////////////////////////////////////////

StgImplementationStd::StgImplementationStd(const std::string& name) :
    StgImplementation(name),
    s_nodes("nodes"),                  // socket sinks
    s_states("states"),                // ...
    s_paststates("pastStates"),        // ...
    s_faceneighcell("faceNeighCell"),  // ...
    s_nvolume("NodalVolume"),          // ...
    s_evolume("ElementVolumes"),       // ...
    s_ienormals("ElementNormals"),     // ...
    s_benormals("BoundariesNormals")   // ...
{
}

//////////////////////////////////////////////////////////////////////////////

void StgImplementationStd::setup()
{
  StgImplementation::setup();
}

//////////////////////////////////////////////////////////////////////////////

void StgImplementationStd::setupVolumes(int& nelm, double& vtot, double& vmin, double& vmax)
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

SafePtr< ConnectivityTable< CFuint > > StgImplementationStd::setupConnectivity()
{
  // return connectivity of first trs found
  std::vector< SafePtr< TopologicalRegionSet > > inner =
    MeshDataStack::getActive()->getFilteredTrsList("inner");
  cf_assert_desc("no trs with inner tag found!",inner.size());
  return inner[0]->getGeo2NodesConn();
}

//////////////////////////////////////////////////////////////////////////////

void StgImplementationStd::setupElementToElement()
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

void StgImplementationStd::PLaS_SetFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp)
{
  CFAUTOTRACE;

  // set number of dimensions
  nbDim = getMethodData().m_plasdata.fp.numDim;

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

void StgImplementationStd::PLaS_SetFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp)
{
  SubSystemStatusStack::getActive()->setDT(fp->dtEul);
  fp->time = (double) SubSystemStatusStack::getActive()->getCurrentTime();  // (double) current time
  fp->iter = (int)    SubSystemStatusStack::getActive()->getNbIter();       // (int) current iteration
  fp->writeOutput = 1;  // (int) flag whether the flow solver writes outputs
}

//////////////////////////////////////////////////////////////////////////////

void StgImplementationStd::PLaS_SetPartitioningData(PLAS_PART_DATA *part)
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

int StgImplementationStd::PLaS_GetElmNode(int elm, int enod)
{
  return (*m_geo2nodes)((CFuint) elm,(CFuint) enod);
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationStd::PLaS_GetElmNeighbour(int elm, int eface)
{
  return m_geo2geo(elm,eface);
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetBndFaceRefCoord(int bnd, int bface, int dim)
{
  return PLaS_GetNodCoord(
    (int) (*getMethodData().m_boundary[bnd]->getGeo2NodesConn())(bface,0),
    dim );
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationStd::PLaS_GetBndDomElm(int bnd, int bface, int dummy)
{
  return (int) s_faceneighcell.getDataHandle()[
    getMethodData().m_boundary[bnd]->getLocalGeoID((CFuint) bface) ].first;
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetNodCoord(int nod, int dim)
{
  return (*(s_nodes.getDataHandle())[nod])[dim];
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetElmNormComp(int elm, int eface, int dim)
{
  return s_ienormals.getDataHandle()[elm][dim][eface];
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetBndFaceNormComp(int bnd, int bface, int dim)
{
  return s_benormals.getDataHandle()[bnd][dim][bface];
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetElmFaceMiddlePoint(int elm, int eface, int dim)
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

double StgImplementationStd::PLaS_GetNodVol(int nod)
{
  return s_nvolume.getDataHandle()[nod];
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetElmVol(int elm)
{
  return s_evolume.getDataHandle()[elm];
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationStd::PLaS_GetNumBndFaces(int bnd)
{
  return (int) getMethodData().m_boundary[bnd]->getGeo2NodesConn()->nbRows();
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationStd::PLaS_GetWallBndFlag(int bnd)
{
  return (getMethodData().m_boundary_is_wall[bnd]? 1:0);
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationStd::PLaS_GetPerBndFlag(int bnd)
{
  return 0;  // no periodic boundaries supported
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetPerBndOffset(int bnd, int dim)
{
  return 0.;  // no periodic boundaries supported
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationStd::PLaS_GetElementType(int elm)
{
  return ELM_SIMPLEX;  // no non-simplex elements supported
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetVelocityComp(int nod, int dim)
{
  if (dim<0 || dim>=nbDim)
    return 0.;
  const int i = getMethodData().m_vidx[dim];
  return (i<0? getMethodData().m_vdef[dim] : (*(s_states.getDataHandle())[nod])[i] );
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetVelocityCompOld(int nod, int dim)
{
  if (dim<0 || dim>=nbDim)
    return 0.;
  const int i = getMethodData().m_vidx[dim];
  return (i<0? getMethodData().m_vdef[dim] : (*(s_paststates.getDataHandle())[nod])[i] );
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetTemperature(int nod)
{
  const int i = getMethodData().m_tidx;
  return (i<0? getMethodData().m_tdef : (*(s_states.getDataHandle())[nod])[i] );
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetTemperatureOld(int nod)
{
  const int i = getMethodData().m_tidx;
  return (i<0? getMethodData().m_tdef : (*(s_paststates.getDataHandle())[nod])[i] );
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetPressure(int nod)
{
  const int i = getMethodData().m_pidx;
  return (i<0? getMethodData().m_pdef : (*(s_states.getDataHandle())[nod])[i] );
}

//////////////////////////////////////////////////////////////////////////////

double StgImplementationStd::PLaS_GetPressureOld(int nod)
{
  const int i = getMethodData().m_pidx;
  return (i<0? getMethodData().m_pdef : (*(s_paststates.getDataHandle())[nod])[i] );
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationStd::PLaS_StartElementSearch(double *pos)
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

int StgImplementationStd::PLaS_EndElementSearch(double *pos)
{
  return (int) m_geo2nodes->nbRows()-1;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

