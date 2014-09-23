
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/SystemFlow.hh"
#include "Muffin/SystemTurb.hh"
#include "Muffin/BCWall.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BCWall,MuffinData,MuffinModule > cBCWall("Wall");

//////////////////////////////////////////////////////////////////////////////

BCWall::BCWall(const std::string& name) :
    BC(name),
    s_mn_volume("NodalVolume"),              // socket sinks
    s_mn_bnormal("NodalNormals"),            // ...
    s_mn_walldistance("NodalWallDistance"),  // ...
    s_mn_wallnode("NodalWallNode")           // ...
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::BCWall");
  addConfigOptionsTo(this);

  // set parameters
  m_walltype = "Static";
  m_temperature = 0.;
  setParameter("Type",&m_walltype);
  setParameter("Temperature",&m_temperature);

  m_zero_vnormal  = true;
  m_zero_vtangent = true;
  m_velocity.clear();
  setParameter("Velocity",&m_velocity);
  setParameter("ZeroVelocityNormal",&m_zero_vnormal);
  setParameter("ZeroVelocityTangent",&m_zero_vtangent);
}

//////////////////////////////////////////////////////////////////////////////

BCWall::~BCWall()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCWall::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;
  options.addConfigOption< std::string >("Type","Wall type, which can be one of \"Static\" (default) or \"StaticTFlux\" (req. temperature present)");
  options.addConfigOption< double >("Temperature","Imposed temperature value or flux (default \"0.\")");

  options.addConfigOption< std::vector< std::string > >("Velocity","Function definition of the forced velocity at the boundaries (default <>)");
  options.addConfigOption< bool >("ZeroVelocityNormal","If normal velocity component should be zero (default true)");
  options.addConfigOption< bool >("ZeroVelocityTangent","If tangencial velocity component should be zero (default true)");
}

//////////////////////////////////////////////////////////////////////////////

void BCWall::setup()
{
  CFAUTOTRACE;
  BC::setup();

  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< RealVector > h_mn_bnormal = s_mn_bnormal.getDataHandle();
  DataHandle< CFreal > h_mn_walldistance = s_mn_walldistance.getDataHandle();

  const CFuint nbDim = (h_nodes.size()?  h_nodes[0]->size()  : 0);
  istflux = (m_walltype=="StaticTFlux");

  RealVector vel(0.,nbDim);          // the velocity to impose
  const RealVector zerov(0.,nbDim);  // auxiliary zero velocity

  if (!m_velocity.size())
    m_velocity.assign(nbDim,"0.");
  if (m_velocity.size()!=nbDim)
    err("\"Velocity\" vector size must be number of dimensions");

  // VectorialFunction setup
  const std::vector< std::string > vf_names(getMethodData().getNodalVariables());
  m_function.setVariables(vf_names);
  RealVector vf_values(vf_names.size());
  m_function.setFunctions(m_velocity);
  try {
    m_function.parse();
  }
  catch (ParserException& e) {
    log("VectorialFunction parsing: " + std::string(e.what()));
    throw;
  }


  // for all the boundary regions, cycle the nodes
  std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
  for (CFuint i=0; i<trs.size(); ++i) {
    const std::vector< CFuint >& nodes = *trs[i]->getNodesInTrs();
    for (std::vector< CFuint >::const_iterator
         n=nodes.begin(); n!=nodes.end(); ++n) {

      // set velocity
      if (!m_zero_vnormal || !m_zero_vtangent) {
        getMethodData().getNodalValues(*n,vf_values);
        m_function.evaluate(vf_values,vel);     // evaluate velocity at node
        RealVector veln(h_mn_bnormal[*n]);      // decompose into
        veln.proj(vel);                         // ... normal and
        RealVector velt(vel-veln);              // ... tangential components
        vel = (m_zero_vtangent? zerov : velt);  // recompose with tangential
        vel += (m_zero_vnormal? zerov : veln);  // ... and normal components
      }
      setInitialState(*n,VVELOCITY,vel);

      // set k and w (Wilcox, p. 202)
      const CFreal wd = (h_mn_walldistance.size()? h_mn_walldistance[*n] : 0.);
      const CFreal w = (wd>0. && Cw2>0.? 6.*nulam/(Cw2*wd*wd) : 0.);  /*FIXME try 60.*/
      setInitialState(*n,VTURBK,1.e-20);
      setInitialState(*n,VTURBW,w);

      // set temperature
      if (!istflux)
        setInitialState(*n,VTEMPERATURE,0.);

    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCWall::applyOnSystemFlow(const SafePtr< System > s, const SafePtr< TopologicalRegionSet > t)
{
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  DataHandle< RealVector > h_mn_bnormal = s_mn_bnormal.getDataHandle();

  const CFuint nbDim = (h_nodes.size()?  h_nodes[0]->size() : 0);
  RealVector vel(0.,nbDim);          // the velocity to impose
  const RealVector zerov(0.,nbDim);  // auxiliary zero velocity
  RealVector vf_values(getMethodData().getNodalVariables().size());


  // for all the boundary regions, cycle the nodes
  const std::vector< CFuint >& nodes = *t->getNodesInTrs();
  for (std::vector< CFuint >::const_iterator
       n=nodes.begin(); n!=nodes.end(); ++n) {

    // set velocity
    if (!m_zero_vnormal || !m_zero_vtangent) {
      getMethodData().getNodalValues(*n,vf_values);
      m_function.evaluate(vf_values,vel);     // evaluate velocity at node
      RealVector veln(h_mn_bnormal[*n]);      // decompose into
      veln.proj(vel);                         // ... normal and
      RealVector velt(vel-veln);              // ... tangential components
      vel = (m_zero_vtangent? zerov : velt);  // recompose with tangential
      vel += (m_zero_vnormal? zerov : veln);  // ... and normal components
    }
    for (CFuint d=0; d<nbDim; ++d)
      setDirichletCondition(*s,*n,1+d, (*h_states[*n])[s->iv+1+d]-vel[d]);

    // set temperature
    if (s.d_castTo< SystemFlow >()->m_coupletemp)
      setDirichletCondition(*s,*n,s->Nsys-1, 0.);

  }
}

//////////////////////////////////////////////////////////////////////////////

void BCWall::applyOnSystemTemp(const SafePtr< System > s, const SafePtr< TopologicalRegionSet > t)
{
  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< RealVector > h_mn_bnormal = s_mn_bnormal.getDataHandle();

  // cycle the nodes
  const std::vector< CFuint >& nodes = *t->getNodesInTrs();
  for (std::vector< CFuint >::const_iterator n=nodes.begin(); n!=nodes.end(); ++n) {

    if (!istflux) {
      setDirichletCondition(*s,*n,0, 0.);
    }
    else {
      h_rhs(*n,s->iv,Neqns) += h_mn_bnormal[*n].sqrNorm() * m_temperature;
    }

  }
}

//////////////////////////////////////////////////////////////////////////////

void BCWall::applyOnSystemTurb(const SafePtr< System > s, const SafePtr< TopologicalRegionSet > t)
{
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  DataHandle< CFreal > h_mn_volume = s_mn_volume.getDataHandle();
  DataHandle< CFreal > h_mn_walldistance = s_mn_walldistance.getDataHandle();
  DataHandle< CFuint > h_mn_wallnode = s_mn_wallnode.getDataHandle();

  // cycle the nodes
  // set k residuals on wall nodes to zero for all models
  const std::vector< CFuint >& nodes = *t->getNodesInTrs();
  for (std::vector< CFuint >::const_iterator n=nodes.begin(); n!=nodes.end(); ++n)
    setDirichletCondition(*s,*n,0, 0.);

  // turbulence model-dependent application (main bc groups)
  switch (turmod) {

    case ITKENA:
    case ITKELB:
    {
      // low re k-e models
      for (std::vector< CFuint >::const_iterator n=nodes.begin(); n!=nodes.end(); ++n) {
        const CFuint intw = h_mn_wallnode[*n];

        const CFreal wd = h_mn_walldistance[*n];
        const CFreal k = (*h_states[intw])[s->iv+0];

        setDirichletCondition( *s,
          *n,1,  (*h_states[*n])[s->iv+1] - 2.*nulam*k/(wd*wd),
          intw,0, -2.*nulam/(wd*wd),
          -h_mn_volume[*n] );

      }
    }
    break;

    case ITKWHR:
    case ITKWLR:
    case ITKWPD:
    case ITKWSS:
    case ITKWBS:
    {
      // high re k-e model and k-w models
      for (std::vector< CFuint >::const_iterator n=nodes.begin(); n!=nodes.end(); ++n)
        setDirichletCondition(*s,*n,1, 0.);
    }
    break;

    case ITKE2L:
    case ITNULL:
    default:
    break;
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

