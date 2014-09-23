
#include "Environment/DirPaths.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/Method.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/Loop.hh"
#include "Muffin/SystemTurb.hh"
#include "Muffin/CC.hh"
#include "Muffin/BC.hh"
#include "Muffin/MuffinData.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

/// Global variables
SafePtr< ConnectivityTable< CFuint > > geo2nodes;  // connectivity list
std::vector< std::vector< int > > neighbors;       // bc nodes neighbors

int Neqns;    // nb. equations in PhysicalModel
int Ndim;     // nb. dimensions in PhysicalModel
int Nvtcell;  // nb. vertices per cell
int Nvtfce;   // nb. vertices per face

turmod_type turmod;
bool turmod_ke;
bool turmod_walldistance;
int buoyancy;
int diffusion;
double nulam;

/// Turbulence model constants
double A1;
double Aeps1;
double Beta;
double Beta1;
double Beta2;
double Betas;
double Bmu;
double Ceps1;
double Ceps2;
double Ck;
double Cmu;
double Cw;
double Cw1;
double Cw2;
double Dmu;
double Gamma1;
double Gamma2;
double Sig1;
double Sig2;
double Sigk1;
double Sigk2;
double Sigw1;
double Sigw2;

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< NullMethodCommand< MuffinData >,MuffinData,MuffinModule > nullMuffinComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void MuffinData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("IntegratorOrder","Order of the Integration to be used for numerical quadrature");
  options.addConfigOption< std::string >("IntegratorQuadrature","Type of Quadrature to be used in the Integration");

  // fluid properties, terms, equations and models
  options.addConfigOption< CFreal >("Density","Fluid density (default 1.0)");
  options.addConfigOption< CFreal >("KinematicViscosity","Kinematic viscosity (default 1.075e-4)");
  options.addConfigOption< CFreal >("SpecificHeatCoeff","Specific heat coefficent [J.g^-1.K^-1] (default 4.1855, water (25 deg C))");
  options.addConfigOption< CFreal >("ThermalConductivity","Thermal conductivity [W.m^-1.K^-1] (default 0.6, water)");
  options.addConfigOption< CFreal >("PrandtlSchmidtNr","Prandtl or Schmidt number (default = 0.7)");
  options.addConfigOption< std::string >("TurbulenceModel","Turbulence model, which can be \"none\" (laminar, default), K-Epsilon \"KEHR\" (high Reynolds - standard wall), K-Epsilon \"KEHG\" (high Reynolds - generalized wall), K-Epsilon \"KE2L\" (two-layer), K-Epsilon \"KELB\" (Lam-Bremhorst), K-Epsilon \"KENA\" (Abe-Kondoh-Nagano), K-Epsilon \"KEV2\" (Durbin K-Epsilon V2), K-Omega \"KWWF\", K-Omega \"KWHR\", K-Omega \"KWLR\", K-Omega \"KWPD\", Menter \"SST\" (SST) or Menter \"BSL\" (baseline)");
  options.addConfigOption< bool >("Diffusion","If diffusion is active (default = true)");
  options.addConfigOption< bool >("Buoyancy","If buoyancy source term is included (default = false)");

  // experimental
  options.addConfigOption< bool >("WallDistance","If wall distance should be calculated (default true)");
  options.addConfigOption< bool >("Restart","Option to restart from the solution provided (default false)");
  options.addConfigOption< std::vector< CFreal > >("Gravity","Gravity vector (default 0. for all dimensions)");
  options.addConfigOption< std::vector< CFreal > >("VelocityReferenceValues","Velocity reference values vector (default < 0. ... >)");
  options.addConfigOption< std::vector< std::string > >("InitialValues","Function definition of the initial values (default < 0. ... >)");
  options.addConfigOption< std::vector< std::string > >("VariableTypes","Variables default types (UNKNOWN, PRESSURE, VELOCITY, SCALAR, TEMPERATURE, CONCENTRATION, POTENTIAL, TURBK, TURBE, TURBW or MAGFIELD) (default < UNKNOWN. ... >)");
}

//////////////////////////////////////////////////////////////////////////////

MuffinData::MuffinData(SafePtr< Method > owner) :
  SpaceMethodData(owner),
  m_convergenceMtd(),
  m_stdTrsGeoBuilder()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_integratorOrderStr = "P1";
  m_integratorQuadratureStr = "INVALID";
  setParameter("IntegratorOrder",&m_integratorOrderStr);
  setParameter("IntegratorQuadrature",&m_integratorQuadratureStr);

  // fluid properties, terms, equations and models
  m_density         = 1.;
  m_kviscosity      = 1.075e-4;
  m_cp              = 4.1855;
  m_k               = 0.6;
  m_prandtl         = 0.7;
  m_turbulencemodel = "none";
  m_diffusion       = true;
  m_buoyancy        = false;
  setParameter("Density",&m_density);
  setParameter("KinematicViscosity",&m_kviscosity);
  setParameter("SpecificHeatCoeff",&m_cp);
  setParameter("ThermalConductivity",&m_k);
  setParameter("PrandtlSchmidtNr",&m_prandtl);
  setParameter("TurbulenceModel",&m_turbulencemodel);
  setParameter("Diffusion",&m_diffusion);
  setParameter("Buoyancy",&m_buoyancy);

  // experimental
  m_wall_distance = false;
  m_restart = false;
  m_gravity.clear();
  m_refvelocity.clear();
  setParameter("WallDistance",&m_wall_distance);
  setParameter("Restart",&m_restart);
  setParameter("Gravity",&m_gravity);
  setParameter("VelocityReferenceValues",&m_refvelocity);

  // initial values function and default variable types definition
  m_initialvalues_def.clear();
  m_vartypes_str.clear();
  setParameter("InitialValues",&m_initialvalues_def);
  setParameter("VariableTypes",&m_vartypes_str);
}

//////////////////////////////////////////////////////////////////////////////

MuffinData::~MuffinData()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void MuffinData::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  SpaceMethodData::configure(args);

  /*
   * setup here the specific integrators for each element with different set
   * of shape and interpolator type
   */
  const CFQuadrature::Type quad  = CFQuadrature::Convert::to_enum(m_integratorQuadratureStr);
  const CFPolyOrder::Type  order = CFPolyOrder::Convert::to_enum(m_integratorOrderStr);
  m_volumeIntegrator.setIntegrationForAllGeo(quad,order);
}

//////////////////////////////////////////////////////////////////////////////

SafePtr< VolumeIntegrator > MuffinData::getVolumeIntegrator()
{
  return &m_volumeIntegrator;
}

//////////////////////////////////////////////////////////////////////////////

void MuffinData::setup()
{
  CFAUTOTRACE;

  // get dimensions and number of equations
  SafePtr< PhysicalModel > pm = PhysicalModelStack::getActive();
  SafePtr< BaseTerm > pm_ct = pm->getImplementor()->getConvectiveTerm();
  Ndim  = pm->getDim();
  Neqns = pm->getNbEq();

  // get variables names (weird stuff!)
  //FIXME hardcoded provider!
  //FIXME isn't this stupid because there will temporarily exist multiple ConvectiveVarSet?
  SelfRegistPtr< ConvectiveVarSet > updatevarset = Environment::Factory< ConvectiveVarSet >::getInstance().getProvider("PhysicalModelDummyPrim")->create(pm_ct);
  cf_assert(updatevarset.isNotNull());
  m_varnames = updatevarset->getVarNames();
  cf_assert(!m_varnames.empty());

  // set default variable types (might be overwritten by the System's)
  if (!m_vartypes_str.size())
    m_vartypes_str.assign(Neqns,"UNKNOWN");
  else if ((int) m_vartypes_str.size()!=Neqns)
    err("\"VariableTypes\" vector size should be number of variables");
  m_vartypes.resize(m_vartypes_str.size());
  for (CFuint i=0; i<m_vartypes.size(); ++i) {
    m_vartypes[i] = (m_vartypes_str[i]=="PRESSURE"?      VPRESSURE :
                    (m_vartypes_str[i]=="VELOCITY"?      VVELOCITY :
                    (m_vartypes_str[i]=="SCALAR"?        VSCALAR :
                    (m_vartypes_str[i]=="TEMPERATURE"?   VTEMPERATURE :
                    (m_vartypes_str[i]=="CONCENTRATION"? VCONCENTRATION :
                    (m_vartypes_str[i]=="POTENTIAL"?     VPOTENTIAL :
                    (m_vartypes_str[i]=="TURBK"?         VTURBK :
                    (m_vartypes_str[i]=="TURBE"?         VTURBE :
                    (m_vartypes_str[i]=="TURBW"?         VTURBW :
                    (m_vartypes_str[i]=="MAGFIELD"?      VMAGFIELD :
                                                         VUNKNOWN ))))))))));
  }

  // set initial values definition vector
  if (!m_initialvalues_def.size())
    m_initialvalues_def.assign(Neqns,"0.");
  else if ((int) m_initialvalues_def.size()!=Neqns && !m_restart)
    err("\"InitialValues\" vector size should be number of variables");

  // set gravity and velocity reference values vector
  if (!m_gravity.size())
    m_gravity.assign(Ndim,0.);
  else if ((int) m_gravity.size()!=Ndim)
    err("\"Gravity\" vector size must be number of dimensions");
  if (!m_refvelocity.size())
    m_refvelocity.assign(Ndim,0.);
  else if ((int) m_refvelocity.size()!=Ndim)
    err("\"VelocityReferenceValues\" vector size must be number of dimensions");

  // set logarithm of L2 error norm for all variables/equations (solution/RHS)
  m_logl2_states.assign(Neqns,0.);
  m_logl2_rhs.assign(Neqns,0.);

  // set SubSystem global status information pointer
  m_status = SubSystemStatusStack::getActive();
}

//////////////////////////////////////////////////////////////////////////////

void MuffinData::setCommands(
  const std::vector< SelfRegistPtr< MuffinCom > >& vcomm_std,
  const std::vector< SelfRegistPtr< MuffinCom > >& vcomm_loops,
  const std::vector< SelfRegistPtr< MuffinCom > >& vcomm_sys,
  const std::vector< SelfRegistPtr< MuffinCom > >& vcomm_cc,
  const std::vector< SelfRegistPtr< MuffinCom > >& vcomm_bc )
{
  std::vector< SelfRegistPtr< MuffinCom > >::const_iterator c;

  m_vcomm_std.clear();
  for (c=vcomm_std.begin(); c!=vcomm_std.end(); ++c)
    m_vcomm_std.push_back(SafePtr< MuffinCom >( c->getPtr() ));

  m_vcomm_loops.clear();
  for (c=vcomm_loops.begin(); c!=vcomm_loops.end(); ++c)
    m_vcomm_loops.push_back(SafePtr< Loop >( c->d_castTo< Loop >().getPtr() ));

  m_vcomm_sys.clear();
  for (c=vcomm_sys.begin(); c!=vcomm_sys.end(); ++c)
    m_vcomm_sys.push_back(SafePtr< System >( c->d_castTo< System >().getPtr() ));

  m_vcomm_cc.clear();
  for (c=vcomm_cc.begin(); c!=vcomm_cc.end(); ++c)
    m_vcomm_cc.push_back(SafePtr< CC >( c->d_castTo< CC >().getPtr() ));

  m_vcomm_bc.clear();
  for (c=vcomm_bc.begin(); c!=vcomm_bc.end(); ++c)
    m_vcomm_bc.push_back(SafePtr< BC >( c->d_castTo< BC >().getPtr() ));
}

//////////////////////////////////////////////////////////////////////////////

void MuffinData::synchronise()
{
  CFAUTOTRACE;
  DataHandle< State*,GLOBAL > h_states = MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();

  PE::GetPE().setBarrier();
  h_states.beginSync();
  h_states.endSync();
  PE::GetPE().setBarrier();
}

//////////////////////////////////////////////////////////////////////////////

int MuffinData::getVariableIndex(const std::string& n)
{
  // return matching variable name index
  cf_assert(!m_varnames.empty());
  for (CFuint i=0; i<m_varnames.size(); ++i) {
    if (m_varnames[i]==n)
      return i;
  }
  ver("variable \"" + n + "\" not found!");
  return -1;
}

//////////////////////////////////////////////////////////////////////////////

int MuffinData::getVariableIndex(const var_type& t)
{
  // return matching variable type index
  cf_assert(!m_vartypes.empty());
  for (CFuint i=0; i<m_vartypes.size(); ++i) {
    if (m_vartypes[i]==t)
      return i;
  }
  ver("variable not found!");
  return -1;
}

//////////////////////////////////////////////////////////////////////////////

const std::string MuffinData::getFilename(const std::string& pre, const std::string& suf, const bool i, const bool t)
{
  // generate filename
  std::ostringstream fn;
  fn << pre;
  if (i)
    fn << "-iter_" << m_status->getNbIter();
  if (t)
    fn << "-time_" << m_status->getCurrentTimeDim();
  fn << suf;

  // prepend ResultsDir to filename
  return (Environment::DirPaths::getInstance().getResultsDir()/fn.str()).string();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< double > MuffinData::getDiffusivity()
{
  std::vector< double > D(Neqns,0.);
  for (CFuint i=0; i<m_vcomm_sys.size(); ++i) {
    const int iv   = m_vcomm_sys[i]->iv;
    const int Nsys = m_vcomm_sys[i]->Nsys;
    std::vector< double > diff = m_vcomm_sys[i]->getDiffusivity();
    for (int e=0; e<Nsys; ++e)
      D[iv+e] = diff[e];
  }
  return D;
}

//////////////////////////////////////////////////////////////////////////////

double MuffinData::getTurbulentViscosity(struct local_node_struct *No_loc, double vol)
{
  /// FIXME this is quite a bad thing
  // get turb_viscosity from first found SystemTurb command
  for (CFuint i=0; i<m_vcomm_sys.size(); ++i) {
    if (m_vcomm_sys[i]->hasTag("Muffin::SystemTurb")) {
      SafePtr< SystemTurb > s = m_vcomm_sys[i].d_castTo< SystemTurb >();
      return s->turb_viscosity(No_loc,vol);
    }
  }
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

void MuffinData::cellgeom(int ic, struct local_node_struct *No_loc, double *vol, int *inc_min)
{
  const ConnectivityTable< CFuint >& geo2nodes = *MeshDataStack::getActive()->getTrs("InnerCells")->getGeo2NodesConn();
  const DataHandle< Node*,GLOBAL > h_nodes = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  for (int i=0; i<Nvtcell; ++i) {
    No_loc[i].node = geo2nodes(ic,i);
    for (int d=0; d<Ndim; ++d)
      No_loc[i].norm[d] = 0.;
  }
  *inc_min = 0;

  if (Ndim==2) {
    // triangles

    for (int i=0; i<Nvtcell; ++i) {
      const int j = (i+1)%Nvtcell;
      const int k = (i+2)%Nvtcell;

      // inwards-facing normals (and vectors k->i)
      No_loc[i].norm[0] = (*h_nodes[No_loc[j].node])[1] - (*h_nodes[No_loc[k].node])[1];
      No_loc[i].norm[1] = (*h_nodes[No_loc[k].node])[0] - (*h_nodes[No_loc[j].node])[0];
      if ( ((*h_nodes[No_loc[i].node])[0] - (*h_nodes[No_loc[k].node])[0]) * No_loc[i].norm[0] +
           ((*h_nodes[No_loc[i].node])[1] - (*h_nodes[No_loc[k].node])[1]) * No_loc[i].norm[1] < 0. )
      err("normal not inward-facing!");

      // set square of magnitude and face with smallest area
      No_loc[i].norm2 =
        No_loc[i].norm[0] * No_loc[i].norm[0] +
        No_loc[i].norm[1] * No_loc[i].norm[1];
      if (No_loc[i].norm2 < No_loc[*inc_min].norm2)
        *inc_min = i;
    }

  }
  else {
    // tetrahedra

    // inward-facing normals
    RealVector v23 = *h_nodes[No_loc[2].node] - *h_nodes[No_loc[3].node];
    RealVector v03 = *h_nodes[No_loc[0].node] - *h_nodes[No_loc[3].node];
    RealVector v13 = *h_nodes[No_loc[1].node] - *h_nodes[No_loc[3].node];
    RealVector v01 = *h_nodes[No_loc[0].node] - *h_nodes[No_loc[1].node];
    RealVector v21 = *h_nodes[No_loc[2].node] - *h_nodes[No_loc[1].node];
    RealVector v31(3);
    v31 = 0.;
    v31 -= v13;

    std::vector< RealVector > normal(4,RealVector(3));
    std::vector< CFreal > sqnormal(4,0.);
    MathTools::MathFunctions::crossProd(v23,v13,normal[0]);
    MathTools::MathFunctions::crossProd(v03,v23,normal[1]);
    MathTools::MathFunctions::crossProd(v01,v31,normal[2]);
    MathTools::MathFunctions::crossProd(v21,v01,normal[3]);
    if (MathTools::MathFunctions::innerProd(v03,normal[0]) < 0.)
      err("normal not inward-facing!");

    // halve magnitude, set square of magnitude and face with smallest area
    for (CFuint i=0; i<4; ++i) {
      normal[i] *= .5;
      sqnormal[i] = normal[i].sqrNorm();
      if (sqnormal[i] < sqnormal[*inc_min])
        *inc_min = i;
    }

    // copy back (this should be removed)
    for (int i=0; i<Nvtcell; ++i) {
      for (int d=0; d<3; ++d)
        No_loc[i].norm[d] = normal[i][d];
      No_loc[i].norm2 = sqnormal[i];
    }

  }

  // volume
  *vol = 0.;
  for (int i=0; i<Nvtcell; ++i)
    for (int d=0; d<Ndim; ++d)
      *vol += (*h_nodes[No_loc[i].node])[d] * No_loc[i].norm[d];
  *vol /= (double) Ndim * (double) Nvtfce;
  if (*vol<0.)
    err("negative volume on element!");
}

//////////////////////////////////////////////////////////////////////////////

void MuffinData::cellgeom(const std::vector< CFuint >& n/*, struct local_node_struct *No_loc, double *vol, int *inc_min*/)
{
  err("time to work");
}

//////////////////////////////////////////////////////////////////////////////

void MuffinData::vecprd3(double *v1, double *v2, double *v)
{
  v[0] = v1[1]*v2[2]-v1[2]*v2[1];
  v[1] = v1[2]*v2[0]-v1[0]*v2[2];
  v[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

//////////////////////////////////////////////////////////////////////////////

void MuffinData::scacde(
  int ic, int iv, struct local_node_struct *No_loc, double vol, int inc_min,
  scalar_convection_type scheme, double diffco, double source, double coeff,
  int coeff_calc )
{
  double u[3];
  double k[4];
  double kplus[4];
  double kminus[4];
  double Win;
  double sumkplus;
  double sumkabs;
  double beta[4];
  double convres;
  double diffres;
  double res[4];
  double q;
  double ninj;

  for (int inc=0; inc<Nvtcell; ++inc)
    No_loc[inc].Res[iv] = 0.;


  // cell-average velocity
  q = 0.;
  for (int id=0; id<Ndim; ++id) {
    u[id] = 0.;
    for (int inc=0; inc<Nvtcell; ++inc)
      u[id] += No_loc[inc].W[id+1];
    u[id] = u[id] / (double) Nvtcell;
    q += u[id]*u[id];
  }
  q = sqrt(q);


  // k[i] = a.dot.n[i]/Nvtfce and convective Dt
  sumkabs  = 0.;
  sumkplus = 0.;
  Win      = 0.;
  for (int inc=0; inc<Nvtcell; ++inc) {
    k[inc] = 0.;
    for (int id=0; id<Ndim; ++id)
      k[inc] += u[id]*No_loc[inc].norm[id];
    k[inc] = k[inc] / (double) Nvtfce;

    kplus[inc] = std::max(k[inc],0.);
    kminus[inc] = std::min(k[inc],0.);

    sumkabs += fabs(k[inc]);
    sumkplus += kplus[inc];

    Win += No_loc[inc].W[iv] * kminus[inc];
  }
  Win /= -sumkplus;


#if 0
  // see if scheme needs changing
  {
    // cell length-scale = Vc/(length or area of smallest edge)
    const double h = 2.*vol*(q<1.e-10?
      1./sqrt(No_loc[inc_min].norm2) :
      q/sumkabs );

    // peclet number
    const double peclet = q*h/diffco;

    if (peclet<20. && peclet>2.)
      scheme = ISSLWS;
    else if (peclet<2.)
      scheme = ISSGAL;
  }
#endif


  // convection

  convres = 0.;
  for (int inc=0; inc<Nvtcell; ++inc)
    convres -= No_loc[inc].W[iv]*k[inc];

  if (scheme==ISSNUL || scheme==ISSFOU || scheme==ISSNSC)
    for (int inc=0; inc<Nvtcell; ++inc)
      No_loc[inc].Res[iv] += source/(double) Nvtcell;
  else
    convres += source;


  // calculate convective nodal residual and coefficents of schemes

  switch (scheme) {

    case ISSFOU:  /* FOU finite-volume scheme */
    {
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int i=1; i<Nvtcell; ++i) {
        const int jnc = (inc+i)%Nvtcell;
        No_loc[inc].Res[iv] += std::max(k[inc]-k[jnc],0.) *
          (No_loc[jnc].W[iv]-No_loc[inc].W[iv])/(double) Nvtcell;
      }
    }
    break;

    case ISSNSC:  /* N-scheme */
    {
      for (int inc=0; inc<Nvtcell; ++inc) {
        res[inc] = kplus[inc] * (No_loc[inc].W[iv] - Win);
        No_loc[inc].Res[iv] -= res[inc];
        if (coeff_calc) {
          No_loc[inc].C[inc] = -kplus[inc];
          for (int i=1; i<Nvtcell; ++i) {
            const int jnc = (inc+i)%Nvtcell;
            No_loc[jnc].C[inc] = -kplus[jnc]*kminus[inc]/sumkplus;
            No_loc[inc].C[jnc] = -kplus[inc]*kminus[jnc]/sumkplus;
          }
        }
      }
    }
    break;

    case ISSGAL:  /* Galerkin scheme */
    {
      for (int inc=0; inc<Nvtcell; ++inc) {
        beta[inc] = 1./(double) Nvtcell;
        No_loc[inc].Res[iv] += beta[inc]*convres;
      }

      if (coeff_calc) {
        for (int inc=0; inc<Nvtcell; ++inc) {
          No_loc[inc].C[inc] = -beta[inc]*k[inc];
          for (int i=1; i<Nvtcell; ++i) {
            const int jnc = (inc+i)%Nvtcell;
            No_loc[jnc].C[inc] = -beta[jnc]*k[inc];
            No_loc[inc].C[jnc] = -beta[inc]*k[jnc];
          }
        }
      }
    }
    break;

    case ISSLDA:  /* LDA-scheme */
    {
      double recipsum = 0.;
      if (sumkplus>PhysicalConstants::_eps)
        recipsum = 1./sumkplus;

      for (int inc=0; inc<Nvtcell; ++inc) {
        beta[inc] = kplus[inc]*recipsum;
        No_loc[inc].Res[iv] += beta[inc]*convres;
      }

      if (coeff_calc) {
        for (int inc=0; inc<Nvtcell; ++inc) {
          No_loc[inc].C[inc] = -beta[inc]*k[inc];
          for (int i=1; i<Nvtcell; ++i) {
            const int jnc = (inc+i)%Nvtcell;
            No_loc[jnc].C[inc] = -beta[jnc]*k[inc];
            No_loc[inc].C[jnc] = -beta[inc]*k[jnc];
          }
        }
      }
    }
    break;

    case ISSLWS:  /* Lax-Wendroff */
    {
      for (int inc=0; inc<Nvtcell; ++inc) {
        beta[inc] = 1./(double) Nvtcell + 0.5*0.5/sumkplus*k[inc];
        No_loc[inc].Res[iv] += beta[inc]*convres;
      }

      if (coeff_calc) {
        for (int inc=0; inc<Nvtcell; ++inc) {
          No_loc[inc].C[inc] = -beta[inc]*k[inc];
          for (int i=1; i<Nvtcell; ++i) {
            const int jnc = (inc+i)%Nvtcell;
            No_loc[jnc].C[inc] = -beta[jnc]*k[inc];
            No_loc[inc].C[jnc] = -beta[inc]*k[jnc];
          }
        }
      }
    }
    break;

    case ISSPSI:  /* PSI-scheme */
    {
      double recipsum = 0.;
      for (int inc=0; inc<Nvtcell; ++inc) {
        res[inc] = No_loc[inc].W[iv] - Win;
        recipsum += kplus[inc]*std::min(0.,res[inc]*convres);
      }
      recipsum = 1./recipsum;

      for (int inc=0; inc<Nvtcell; ++inc) {
        beta[inc] = kplus[inc]*std::min(0.,res[inc]*convres)*recipsum;
        No_loc[inc].Res[iv] += beta[inc]*convres;
      }

      if (coeff_calc) {
        for (int inc=0; inc<Nvtcell; ++inc) {
          /*
           * C[inc][inc] = -beta[inc]*k[inc];
           * for (int i=1; i<Nvtcell; ++i) {
           *   const int jnc = (inc+i)%Nvtcell;
           *   C[jnc][inc] = -beta[jnc]*k[inc];
           *   C[inc][jnc] = -beta[inc]*k[jnc];
           * }
           */
          /* N-scheme Jacobian entries */
          No_loc[inc].C[inc] = -kplus[inc];
          for (int i=1; i<Nvtcell; ++i) {
            const int jnc = (inc+i)%Nvtcell;
            No_loc[jnc].C[inc] = -kplus[jnc]*kminus[inc]/sumkplus;
            No_loc[inc].C[jnc] = -kplus[inc]*kminus[jnc]/sumkplus;
          }
        }
      }

      for (int inc=0; inc<Nvtcell; ++inc)
        beta[inc] = kplus[inc]/sumkplus;
    }
    break;

    case ISSNUL:  /* No convection */
    default:
    {
      for (int inc=0; inc<Nvtcell; ++inc)
        for (int jnc=0; jnc<Nvtcell; ++jnc)
          No_loc[inc].C[jnc] = 0.;
    }
    break;

  }


  /* diffusive residual */
  if (diffusion) {

    // gradient of variable
    RealVector grad(Ndim);
    grad = 0.;
    for (int id=0; id<Ndim; ++id)
      for (int inc=0; inc<Nvtcell; ++inc)
        grad[id] += No_loc[inc].W[iv] * No_loc[inc].norm[id];
    grad /= (double) Nvtfce*vol;

    // distribution
    for (int inc=0; inc<Nvtcell; ++inc) {
      diffres = 0.;
      for (int id=0; id<Ndim; ++id)
        diffres += grad[id]*No_loc[inc].norm[id];
      diffres *= -diffco/(double) Nvtfce;

      No_loc[inc].Res[iv] += diffres;
    }
  }


  if (coeff_calc) {

    /* add source terms */
    if (scheme==ISSNUL || scheme==ISSFOU || scheme==ISSNSC) {
      for (int inc=0; inc<Nvtcell; ++inc)
        for (int jnc=0; jnc<Nvtcell; ++jnc)
          No_loc[inc].C[jnc] += coeff/(double) Nvtcell;
    }
    else {
      for (int inc=0; inc<Nvtcell; ++inc)
        for (int jnc=0; jnc<Nvtcell; ++jnc)
          No_loc[inc].C[jnc] += coeff*beta[inc];
    }

    /* add diffusion terms */
    if (diffusion) {
      const double diffac = diffco/((double) Nvtfce*(double) Nvtfce*vol);
      for (int inc=0; inc<Nvtcell; ++inc) {
        for (int jnc=0; jnc<Nvtcell; ++jnc) {
          ninj = 0.;
          for (int id=0; id<Ndim; ++id)
            ninj += No_loc[inc].norm[id]*No_loc[jnc].norm[id];
          No_loc[inc].C[jnc] -= diffac*ninj;
        }
      }
    }

  }
}

//////////////////////////////////////////////////////////////////////////////

bool MuffinData::ispointintriangle(
  const double x,  const double y,
  const double x1, const double y1,
  const double x2, const double y2,
  const double x3, const double y3 )
{
  // vectors
  double v0[2] = { x3-x1, y3-y1 };
  double v1[2] = { x2-x1, y2-y1 };
  double v2[2] = { x -x1, y -y1 };

  // dot products
  const double dot00 = v0[0]*v0[0] + v0[1]*v0[1];
  const double dot01 = v0[0]*v1[0] + v0[1]*v1[1];
  const double dot02 = v0[0]*v2[0] + v0[1]*v2[1];
  const double dot11 = v1[0]*v1[0] + v1[1]*v1[1];
  const double dot12 = v1[0]*v2[0] + v1[1]*v2[1];

  // barycentric coordinates
  const double invdenom = 1. / (dot00*dot11 - dot01*dot01);
  const double u = (dot11*dot02 - dot01*dot12)*invdenom;
  const double v = (dot00*dot12 - dot01*dot02)*invdenom;

  // check if point is in triangle
  return ((u>=0) && (v>=0) && (u+v<1+PhysicalConstants::_eps));
}

//////////////////////////////////////////////////////////////////////////////

double MuffinData::interpolate(
  const double x,  const double y,
  const double x1, const double y1, const double v1,
  const double x2, const double y2, const double v2,
  const double x3, const double y3, const double v3 )
{
  const double nx =  (y2-y1)*(v3-v1)-(y3-y1)*(v2-v1);
  const double ny = -(x2-x1)*(v3-v1)+(x3-x1)*(v2-v1);
  double nz = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
  nz = (nz>0?
    std::max(nz,PhysicalConstants::_eps) :
    std::min(nz,-PhysicalConstants::_eps) );
  return (nx*(x1-x) + ny*(y1-y))/nz + v1;
}

//////////////////////////////////////////////////////////////////////////////

double MuffinData::distance(
  const double x1, const double y1, const double z1,
  const double x2, const double y2, const double z2 )
{
  return sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) );
}

//////////////////////////////////////////////////////////////////////////////

double MuffinData::distanceFaceLine2D(
  const double x,  const double y,
  const double x1, const double y1,
  const double x2, const double y2 )
{
  const double dx = x2-x1;
  const double dy = y2-y1;
  return fabs( dx*(y1-y)-dy*(x1-x) )/sqrt( dx*dx + dy*dy );
}

//////////////////////////////////////////////////////////////////////////////

double MuffinData::distanceFacePlane3D(
  const double x,  const double y,  const double z,
  const double x1, const double y1, const double z1,
  const double x2, const double y2, const double z2,
  const double x3, const double y3, const double z3 )
{
  const double nx =  (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
  const double ny = -(x2-x1)*(z3-z1)-(x3-x1)*(z2-z1);
  const double nz =  (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
  const double nn = sqrt( nx*nx + ny*ny + nz*nz );
  return fabs( nx*(x-x1) + ny*(y-y1) + nz*(z-z1) )/nn;
}

//////////////////////////////////////////////////////////////////////////////

double** MuffinData::allocate_double_matrix(int nrh, int nch)
{
  /* allocate pointers to rows */
  double **m = new double*[ nrh ];
  if (!m)
    err("allocation failure 1 in allocate_double_matrix");

  /* allocate rows and set pointers to them */
  m[0] = new double[ nrh*nch ];
  if (!m[0])
    err("allocation failure 1 in allocate_double_matrix");
  for (int i=1; i<nrh; ++i)
    m[i] = m[i-1]+nch;

  /* initialize */
  for (int i=0; i<nrh; ++i)
    for (int j=0; j<nch; ++j)
      m[i][j] = 0.;

  /* return pointer to array of pointers to rows */
  return m;
}

//////////////////////////////////////////////////////////////////////////////

double*** MuffinData::allocate_double_tensor(int nrh, int nch, int ndh)
{
  double ***t = new double**[ nrh ];
  for (int i=0; i<nrh; ++i)
    t[i] = allocate_double_matrix(nch,ndh);
  return t;
}

//////////////////////////////////////////////////////////////////////////////

int** MuffinData::allocate_int_matrix(int nrh, int nch)
{
  /* allocate pointers to rows */
  int **m = new int*[ nrh ];
  if (!m)
    err("allocation failure 1 in allocate_int_matrix");

  /* allocate rows and set pointers to them */
  m[0] = new int[ nrh*nch ];
  if (!m[0])
    err("allocation failure 1 in allocate_int_matrix");
  for (int i=1; i<nrh; ++i)
    m[i] = m[i-1]+nch;

  /* initialize */
  for (int i=0; i<nrh; ++i)
    for (int j=0; j<nch; ++j)
      m[i][j] = 0;

  /* return pointer to array of pointers to rows */
  return m;
}

//////////////////////////////////////////////////////////////////////////////

double* MuffinData::dvector(long nl, long nh)
{
  double *v = new double[ nh-nl+1 ];
  if (!v)
    err("allocation failure in dvector");
  return v-nl;
}

//////////////////////////////////////////////////////////////////////////////

double** MuffinData::dmatrix(long nrl, long nrh, long ncl, long nch)
{
  const long nrow = nrh-nrl+1;
  const long ncol = nch-ncl+1;

  /* allocate pointers to rows */
  double **m = new double*[ nrow ];
  if (!m)
    err("allocation failure in dmatrix (rows)");
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = new double[ nrow*ncol ];
  if (!m[nrl])
    err("allocation failure in dmatrix (columns)");
  m[nrl] -= ncl;

  for (long i=nrl+1; i<=nrh; ++i)
    m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

//////////////////////////////////////////////////////////////////////////////

double*** MuffinData::d3tensor(
  long nrl, long nrh, long ncl, long nch, long ndl, long ndh )
{
  const long nrow = nrh-nrl+1;
  const long ncol = nch-ncl+1;
  const long ndep = ndh-ndl+1;

  /* allocate pointers to pointers to rows */
  double ***t = new double**[ nrow ];
  if (!t)
    err("allocation failure in d3tensor (rows)");
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl] = new double*[ nrow*ncol ];
  if (!t[nrl])
    err("allocation failure in d3tensor (columns)");
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl] = new double[ nrow*ncol*ndep ];
  if (!t[nrl][ncl])
    err("allocation failure in d3tensor (depth)");
  t[nrl][ncl] -= ndl;

  for (long j=ncl+1; j<=nch; ++j)
    t[nrl][j]=t[nrl][j-1]+ndep;

  for (long i=nrl+1; i<=nrh; ++i) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for (long j=ncl+1; j<=nch; ++j)
      t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

//////////////////////////////////////////////////////////////////////////////

void MuffinData::free_dvector(double *v, long nl, long nh)
{
  double *xv = v+nl;
  delete[] xv;
}

//////////////////////////////////////////////////////////////////////////////

void MuffinData::free_dmatrix(
  double **m, long nrl, long nrh, long ncl, long nch )
{
  double *xv = m[nrl]+ncl;
  double **xm = m+nrl;
  delete[] xv;
  delete[] xm;
}

//////////////////////////////////////////////////////////////////////////////

void MuffinData::free_d3tensor(
  double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh )
{
  double *xv = t[nrl][ncl]+ndl;
  double **xm = t[nrl]+ncl;
  double ***xt = t+nrl;
  delete[] xv;
  delete[] xm;
  delete[] xt;
}

//////////////////////////////////////////////////////////////////////////////

void MuffinData::getNodalValues(const CFuint n, RealVector& v)
{
  const DataHandle< Node*,GLOBAL > h_nodes = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  const CFuint nbDim = (h_nodes.size()? h_nodes[0]->size():0);

  cf_assert_desc("getNodalValues: unexpected nodal values vector size",v.size()==nbDim+1+1);
  cf_assert_desc("getNodalValues: unexpected node number",n<h_nodes.size());

  // get coordinates, iteration and time values
  for (CFuint d=0; d<nbDim; ++d)
    v[d] = (*h_nodes[n])[d];
  v[nbDim+0] = m_status->getNbIter();
  v[nbDim+1] = m_status->getCurrentTimeDim();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< std::string > MuffinData::getNodalVariables()
{
  const DataHandle< Node*,GLOBAL > h_nodes = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  const CFuint nbDim = (h_nodes.size()? h_nodes[0]->size():0);

  // get coordinates, iteration and time names
  std::vector< std::string > v(nbDim+2);
  for (CFuint d=0; d<nbDim; ++d)
    v[d] = std::string(1,char('X'+d));
  v[nbDim+0] = "i";
  v[nbDim+1] = "t";
  return v;
}

//////////////////////////////////////////////////////////////////////////////

  }  // end of namespace Muffin
}  // end of namespace COOLFluiD

