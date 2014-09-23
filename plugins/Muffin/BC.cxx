
#include "Common/GlobalReduce.hh"
#include "Muffin/System.hh"
#include "Muffin/BC.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

BC::BC(const std::string& name) :
    MuffinCom(name),
    s_rhs("rhs"),                      // socket sinks
    s_nodes("nodes"),                  // ...
    s_states("states"),                // ...
    s_faceneighcell("faceNeighCell"),  // ...
    s_mn_priority("NodalPriority"),    // ...
    m_point_index(-1)
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  addConfigOptionsTo(this);

  m_point_xyz.clear();
  m_bnpriority = 1;
  setParameter("applyPoint",&m_point_xyz);
  setParameter("Priority",&m_bnpriority);
}

//////////////////////////////////////////////////////////////////////////////

BC::~BC()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BC::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector< double > >("applyPoint"," (default <>)");
  options.addConfigOption< CFuint >("Priority","Priority in applicable nodes over other boundary conditions (default 1)");
}

//////////////////////////////////////////////////////////////////////////////

void BC::configure(Config::ConfigArgs& args)
{
  MuffinCom::configure(args);
  if (m_bnpriority<1)
    err("incorrect \"Priority\" option, should be greater than 0!");
}

//////////////////////////////////////////////////////////////////////////////

void BC::setup()
{
  CFAUTOTRACE;
  MuffinCom::setup();

  // setup m_point_index to closest point to m_point_xyz
  if (m_point_xyz.size()) {

    if (m_point_xyz.size()!=(unsigned) Ndim)
      err("\"applyPoint\" needs to be same size as dimensions");

    // seek closest node
    DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
    double distance2_min = 1.e99;
    for (CFuint n=0; n<h_nodes.size(); ++n) {
      if (h_nodes[n]->isParUpdatable()) {
        const double dx = (*h_nodes[n])[0] - m_point_xyz[0];
        const double dy = (*h_nodes[n])[1] - m_point_xyz[1];
        const double dz = (Ndim>2? (*h_nodes[n])[2] - m_point_xyz[2] : 0.);
        const double d2 = dx*dx+dy*dy+dz*dz;
        if (d2<distance2_min) {
          distance2_min = d2;
          m_point_index = (int) n;
        }
      }
    }

    if (m_point_index>=0) {
      log( "bulk application point (i:x,y,z): " + StringOps::to_str(m_point_index) +
        ":" + StringOps::to_str((*h_nodes[m_point_index])[0]) +
        "," + StringOps::to_str((*h_nodes[m_point_index])[1]) +
        "," + StringOps::to_str(Ndim>2? (*h_nodes[m_point_index])[2] : 0.) );

      // get mesh connectivity
        SafePtr< TopologicalRegionSet > trs =
          MeshDataStack::getActive()->getTrs("InnerCells");
        ConnectivityTable< CFuint > geo2nodes = *trs->getGeo2NodesConn();

      // add neighbors (includes itself) to nnodes, removing duplicates
      const CFuint n = (CFuint) m_point_index;
      std::vector< int > nnodes;
      for (CFuint ic=0; ic<trs->getLocalNbGeoEnts(); ++ic)
        if ( geo2nodes(ic,0)==n ||
             geo2nodes(ic,1)==n ||
             geo2nodes(ic,2)==n ||
             (Ndim>2? geo2nodes(ic,3)==n : false ) ) {
          nnodes.push_back(geo2nodes(ic,0));
          nnodes.push_back(geo2nodes(ic,1));
          nnodes.push_back(geo2nodes(ic,2));
          if (Ndim>2)
            nnodes.push_back(geo2nodes(ic,3));
        }
      sort(nnodes.begin(),nnodes.end());
      nnodes.erase( unique(nnodes.begin(),nnodes.end()), nnodes.end() );

      // set neighbors
      neighbors[n] = nnodes;
    }

  }
}

//////////////////////////////////////////////////////////////////////////////

void BC::apply(const SafePtr< System > s)
{
  const std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
  for (std::vector< SafePtr< TopologicalRegionSet > >::const_iterator t=trs.begin(); t!=trs.end(); ++t) {
    if      (s->hasTag("Muffin::SystemFlow"))   { applyOnSystemFlow(s,*t);   }
    else if (s->hasTag("Muffin::SystemTemp"))   { applyOnSystemTemp(s,*t);   }
    else if (s->hasTag("Muffin::SystemTurb"))   { applyOnSystemTurb(s,*t);   }
    else if (s->hasTag("Muffin::SystemMITReM")) { applyOnSystemMITReM(s,*t); }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BC::setInitialState(const CFuint& n, const var_type& state, const CFreal& v)
{
  DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();
  if (h_mn_priority[n]>m_bnpriority)
    return;

  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  for (CFuint i=0; i<getMethodData().m_vartypes.size(); ++i) {
    if (getMethodData().m_vartypes[i]==state)
      (*h_states[n])[i] = v;
  }
}

//////////////////////////////////////////////////////////////////////////////

void BC::setInitialState(const CFuint& n, const var_type& state, const RealVector& vv)
{
  DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();
  if (h_mn_priority[n]>m_bnpriority || !vv.size())
    return;

  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  const std::vector< var_type >& vartypes = getMethodData().m_vartypes;
  for (CFuint i=0; i<vartypes.size(); ++i) {
    if (vartypes[i]==state) {
      for (CFuint d=0; d<vv.size(); ++d)
        (*h_states[n])[i+d] = vv[d];
      i += vv.size()-1;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BC::setDirichletCondition( System& s,
  const CFuint n, const CFuint e, const double v, const double f )
{
  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();
  if (h_mn_priority[n]>m_bnpriority)
    return;

  std::vector< CFint > matrixid(neighbors[n].size(),0);
  CFint matrixid_n = -1;
  const LSSIdxMapping& idxMapping =
    getMethodData().getLinearSystemSolver()[s.is]->getLocalToGlobalMapping();
  for (int nb=0; nb<(int) neighbors[n].size(); ++nb) {
    matrixid[nb] = (CFint) idxMapping.getColID(neighbors[n][nb]);
    if ((CFuint) neighbors[n][nb]==n)
      matrixid_n = matrixid[nb];
  }
  cf_assert_desc("node index not found!",matrixid_n>=0);

  const CFint mrow = matrixid_n*s.Nsys+(CFint) e;
  h_rhs(n,s.iv+e,Neqns) = v*f;
  for (int nb=0; nb<(int) neighbors[n].size(); ++nb) {
    const CFint mcol0 = matrixid[nb]*s.Nsys;
    for (CFint j=0; j<(CFint) s.Nsys; ++j)
      s.matrix->setValue(mrow,mcol0+j, 0. );
  }
  s.matrix->setValue(mrow,mrow, f);
}

//////////////////////////////////////////////////////////////////////////////

void BC::setDirichletCondition( System& s,
  const CFuint n1, const CFuint e1, const double v1,
  const CFuint n2, const CFuint e2, const double v2, const double f )
{
  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();

  std::vector< CFint > matrixid(neighbors[n1].size(),0);
  CFint matrixid_n1 = -1;
  CFint matrixid_n2 = -1;
  const LSSIdxMapping& idxMapping = getMethodData().getLinearSystemSolver()[s.is]->getLocalToGlobalMapping();
  for (int nb=0; nb<(int) neighbors[n1].size(); ++nb) {
    matrixid[nb] = (CFint) idxMapping.getColID(neighbors[n1][nb]);
    if ((CFuint) neighbors[n1][nb]==n1)
      matrixid_n1 = matrixid[nb];
    if ((CFuint) neighbors[n1][nb]==n2)
      matrixid_n2 = matrixid[nb];
  }
  cf_assert_desc("node index not found!",matrixid_n1>=0 && matrixid_n2>=0);

   const CFint mrow = matrixid_n1*s.Nsys+(int) e1;
   h_rhs(n1,s.iv+e1,Neqns) = v1*f;
   for (int nb=0; nb<(int) neighbors[n1].size(); ++nb) {
     const CFint mcol0 = matrixid[nb]*s.Nsys;
     for (CFint j=0; j<(CFint) s.Nsys; ++j)
       s.matrix->setValue(mrow,mcol0+j, 0. );
   }
   const CFint mcol = matrixid_n2*s.Nsys+(int) e2;
   s.matrix->setValue(mrow,mcol, v2*f);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

