
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/System.hh"
#include "Muffin/BCMassFlux.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BCMassFlux,MuffinData,MuffinModule > cBCMassFluxProvider("MassFlux");

//////////////////////////////////////////////////////////////////////////////

BCMassFlux::BCMassFlux(const std::string& name) :
    BC(name),
    m_velocity_i(-1)
{
  attachTag("Muffin::BCMassFlux");
}

//////////////////////////////////////////////////////////////////////////////

void BCMassFlux::setup()
{
  CFAUTOTRACE;
  BC::setup();
  MuffinData& d = getMethodData();


  log("get velocity field...");
  for (CFuint i=0; i<d.m_vartypes.size(); ++i)
    if (d.m_vartypes[i]==VVELOCITY) {
      m_velocity_i = i;
      break;
    }
  ver("get velocity field.");


#if 0
  log("setup mass flux diffusivity coefficients and variables...");
  m_massflux_d.clear();
  m_massflux_i.clear();
  m_massflux_n.clear();
  if (m_velocity_vars.size()) {
    std::vector< double > diff = d.getDiffusivity();
    for (CFuint i=0; i<diff.size(); ++i)
      if (diff[i]>0.) {
        m_massflux_d.push_back(diff[i]);
        m_massflux_i.push_back(i);
        m_massflux_n.push_back(d.m_varnames[i]);
      }
  }
  ver("setup mass flux diffusivity coefficients and variables.");
#endif
}

//////////////////////////////////////////////////////////////////////////////

#if 0
void BCMassFlux::applyOnSystemFlow(System& s)
{
  CFAUTOTRACE;
  MuffinData& data = getMethodData();

  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  DataHandle< std::pair< CFuint,CFuint > > h_faceneighcell = s_faceneighcell.getDataHandle();

  std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
  for (CFuint i=0; i<m_massflux_i.size(); ++i) {
    for (CFuint t=0; t<trs.size(); ++t) {


      double trs_area   = 0.;
      double trs_flux_a = 0.;
      double trs_flux_d = 0.;

      std::vector< double > norm,
                            a, b, c, d,
                            C;
    
      SafePtr< ConnectivityTable< CFuint > > faces = trs[t]->getGeo2NodesConn();
      for (CFuint ifc=0; ifc<faces->nbRows(); ++ifc) {
    
        // inner cell and face index, converted from coolfluid to "opposite node"
        const CFuint _c = h_faceneighcell[ trs[t]->getLocalGeoID(ifc) ].first;
        const CFuint _f = h_faceneighcell[ trs[t]->getLocalGeoID(ifc) ].second;
        const int of = ( Ndim==2?  (_f==0? 2 :
                                    _f==1? 0 :
                                           1 ) :
                                  ((_f==0? 3 :
                                   (_f==1? 2 :
                                    _f==2? 0 :
                                           1 ))) );
    
        // geometrical properties (face indices, plus node oposite to face)
        double area;
        std::vector< CFuint > n(Nvtcell,0);
        for (CFuint v=0; v<Nvtcell-1; ++v)
          n[v] = faces(ifc,v);
        n[Nvtcell-1] = (*geo2nodes)(_c,of);
        if ( !h_nodes[n[0]]->isParUpdatable() ||
             !h_nodes[n[1]]->isParUpdatable() ||
             !h_nodes[n[2]]->isParUpdatable() ||
             Ndim==2? false : !h_nodes[n[3]]->isParUpdatable() )
          continue;
        getGeometry(n,norm,&area,a,b,c,d);
    
        // get velocity field
        std::vector< double > vx(Nvtcell,0.),
                              vy(Nvtcell,0.),
                              vz(Nvtcell,0.);
        for (unsigned i=0; i<Nvtcell; ++i) {
          vx[i] = (*h_states[n[i]])[m_velocity_i+0];
          vy[i] = (*h_states[n[i]])[m_velocity_i+1];
          vz[i] = (*h_states[n[i]])[m_velocity_i+2];
        }
    
        // get concentrations
        std::vector< double > C(Nvtcell,data.m_density);
        if (m_massflux_i[i]>=0) {
          for (unsigned j=0; j<Nvtcell; ++j)
            C[j] = (*h_states[n[j]])[m_massflux_i[i]]
        }
    
        // boundary element advective flux (integral of density in each direction)
        *trs_flux_a += area / 12. * (
            norm[0] * (
            C[0]*vx[0]*2. + C[0]*vx[1]    + C[0]*vx[2] +
            C[1]*vx[0]    + C[1]*vx[1]*2. + C[1]*vx[2] +
            C[2]*vx[0]    + C[2]*vx[1]    + C[2]*vx[2]*2. ) +
            norm[1] * (
            C[0]*vy[0]*2. + C[0]*vy[1]    + C[0]*vy[2] +
            C[1]*vy[0]    + C[1]*vy[1]*2. + C[1]*vy[2] +
            C[2]*vy[0]    + C[2]*vy[1]    + C[2]*vy[2]*2. ) +
            norm[2] * (
            C[0]*vz[0]*2. + C[0]*vz[1]    + C[0]*vz[2] +
            C[1]*vz[0]    + C[1]*vz[1]*2. + C[1]*vz[2] +
            C[2]*vz[0]    + C[2]*vz[1]    + C[2]*vz[2]*2. ) );
    
        // boundary element diffusive flux
        *trs_flux_d += area * m_massflux_d[i] * (
            norm[0] * ( C[0]*b[0] + C[1]*b[1] + C[2]*b[2] + C[3]*b[3] ) +
            norm[1] * ( C[0]*c[0] + C[1]*c[1] + C[2]*c[2] + C[3]*c[3] ) +
            norm[2] * ( C[0]*d[0] + C[1]*d[1] + C[2]*d[2] + C[3]*d[3] ) );
    
        // boundary element area
        *trs_area += area;
      }
    
      // synchronize flux values
      Common::GlobalReduceOperation< Common::GRO_SUM >(trs_flux_a,trs_flux_a);
      Common::GlobalReduceOperation< Common::GRO_SUM >(trs_flux_d,trs_flux_d);

      // log result
      std::string s( "mass flux of \"" + m_massflux_n[i] + "\", TRS \"" +
        trs[t]->getName() + "\" (area: " + StringOps::to_str(trs_area) + "):" );
      s += "  advective: "  + StringOps::to_str(trs_flux_a);
      s += "  diffusive: "  + StringOps::to_str(trs_flux_d);
      s += "  total: " + StringOps::to_str(trs_flux_a+trs_flux_d);
      log(s);


    }
  }
}
#endif

//////////////////////////////////////////////////////////////////////////////

void BCMassFlux::getGeometry(
  std::vector< int >& n, std::vector< double >& norm, double* area,
  std::vector< double >& a, std::vector< double >& b,
  std::vector< double >& c, std::vector< double >& d )
{
  // coordinates
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  const int N = n.size();
  std::vector< double > x(N,0.);
  std::vector< double > y(N,0.);
  std::vector< double > z(N,0.);
  for (int i=0; i<N; ++i) {
    x[i] = (*h_nodes[n[i]])[0];
    y[i] = (*h_nodes[n[i]])[1];
    z[i] = (Ndim>2? (*h_nodes[n[i]])[2] : 0.);
  }

  // normal components, area and volume
  const double nx = (y[2]-y[0]) * (z[1]-z[0]) - (z[2]-z[0]) * (y[1]-y[0]);
  const double ny = (z[2]-z[0]) * (x[1]-x[0]) - (x[2]-x[0]) * (z[1]-z[0]);
  const double nz = (x[2]-x[0]) * (y[1]-y[0]) - (y[2]-y[0]) * (x[1]-x[0]);
  const double _a = sqrt(nx*nx + ny*ny + nz*nz)/2.;
  const double _v = (
    (x[1]-x[0])*( (y[2]-y[0])*(z[3]-z[0]) - (y[3]-y[0])*(z[2]-z[0]) ) -
    (x[2]-x[0])*( (y[1]-y[0])*(z[3]-z[0]) - (y[3]-y[0])*(z[1]-z[0]) ) +
    (x[3]-x[0])*( (y[1]-y[0])*(z[2]-z[0]) - (y[2]-y[0])*(z[1]-z[0]) ) )/6.;
  *area = _a;

  // shape function coefficients
  a.resize(N);
  b.resize(N);
  c.resize(N);
  d.resize(N);
  const int ij[] = {1,0,0,0};
  const int ik[] = {3,2,3,1};
  const int il[] = {2,3,1,2};
  const double iv6 = 1./(6.*_v);
  for (int i=0; i<N; ++i) {
    const int j = ij[i];
    const int k = ik[i];
    const int l = il[i];
    a[i] = 0. /*- iv6 * (
      x[j]*(y[k]*z[l]-y[l]*z[k]) +
      x[k]*(y[l]*z[j]-y[j]*z[l]) +
      x[l]*(y[j]*z[k]-y[k]*z[j]) ) */ ;
    b[i] = iv6*( y[j]*(z[k]-z[l]) + y[k]*(z[l]-z[j]) + y[l]*(z[j]-z[k]) );
    c[i] = iv6*( x[j]*(z[l]-z[k]) + x[k]*(z[j]-z[l]) + x[l]*(z[k]-z[j]) );
    d[i] = iv6*( x[j]*(y[k]-y[l]) + x[k]*(y[l]-y[j]) + x[l]*(y[j]-y[k]) );
  }

  // normalize norm components and set face area
  norm.resize(3);
  norm[0] = nx/(2.*_a);
  norm[1] = ny/(2.*_a);
  norm[2] = nz/(2.*_a);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

