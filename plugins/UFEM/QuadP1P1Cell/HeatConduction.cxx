#include "Environment/ObjectProvider.hh"

#include "UFEM/QuadP1P1Cell/HeatConduction.hh"
#include "UFEM/QuadP1P1Cell/UFEMQuadP1P1Cell.hh"
#include "UFEM/functors/GradSquared.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace QuadP1P1Cell {

using namespace std;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider < HeatConduction,
                              UFEMTerm,
                              UFEMQuadP1P1CellPlugin,
                              UFEMTerm::NARGS >
aQuadP1P1Cell_HeatConduction_Provider ( "QuadP1P1Cell_HeatConduction" );

//////////////////////////////////////////////////////////////////////////////

namespace detail
{

/// Gauss quadrature of a function of xi and eta
template<typename FunctorT, typename ValueT>
void quad(FunctorT& functor, ValueT& result) {
  //static const CFreal gp = sqrt(1./3.); // Gauss points coordinate for 2-point integral
  result = 4. * functor(0., 0.); // one-point integration
  //result = functor(-gp, -gp) + functor(-gp, gp) + functor(gp, -gp) + functor(gp, gp);
}

} // namespace detail

//////////////////////////////////////////////////////////////////////////////

void HeatConduction::defineConfigOptions ( Config::OptionList& options )
{
  CFAUTOTRACE;

  options.addConfigOption< CFreal >("a", "Heat Conduction Coefficient");
  options.addConfigOption< bool >("CFIntegrator","Use the framework numeric integration method");
}

//////////////////////////////////////////////////////////////////////////////

HeatConduction::HeatConduction ( const std::string& name, Common::SafePtr<ElemProps> props  ) :
  UFEMTerm ( name, props ),
  socket_interStates("interStates")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  //Default Effective Heat Conductivity
  m_a     = 1.;

  m_cf_integrator = false;

  setParameter( "a",         &m_a );
  setParameter( "CFIntegrator",         &m_cf_integrator );
}

//////////////////////////////////////////////////////////////////////////////

HeatConduction::~HeatConduction() 
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void HeatConduction::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  UFEMTerm::configure ( args );

  CFLog(INFO, getClassName() << ": heat conductivity: "  << m_a      << "\n");
  CFLog(INFO, getClassName() << ": CFIntegrator: "  << m_cf_integrator << "\n");

}

//////////////////////////////////////////////////////////////////////////////

void HeatConduction::setup ()
{
  CFAUTOTRACE;
  UFEMTerm::setup ();

  socket_interStates.setParentNamespace ( getMethodData().getNamespace() );
}

//////////////////////////////////////////////////////////////////////////////

void HeatConduction::unsetup ()
{
  CFAUTOTRACE;
  UFEMTerm::unsetup ();
}

//////////////////////////////////////////////////////////////////////////////

void HeatConduction::compute ( const Framework::GeometricEntity& cell , AssemblyData& adata )
{
  CFAUTOTRACE;

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nb_states = cell.nbStates();

  CFuint Ti,Tj;
  COOLFluiD::RealMatrix K(nb_states, nb_states);
  VolumeIntegrator integrator;
  integrator.setup();
  functors::GradSquared functor(cell);

  if(m_cf_integrator) {
    integrator.integrateGeometricFunctorOnGeoEnt(const_cast<Framework::GeometricEntity*>(&cell), functor, K);
  } else {
    detail::quad(functor, K);
  }

  K *= m_a;

  Ti=0;
  for(CFuint i=0; i<nb_states; ++i) {

    Tj=0;
    for(CFuint j=0; j<nb_states; ++j) {

      //diffusion (Standard)
      adata.A(Ti,Tj) += K[j + i*nb_states];

      Tj+=nbEqs;
    }
    Ti+=nbEqs;
  }

  return;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > > HeatConduction::providesSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > > HeatConduction::needsSockets()
{
  CFAUTOTRACE;
  std::vector< SafePtr< BaseDataSocketSink > > result = UFEMTerm::needsSockets();

  result.push_back(&socket_interStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    }  // namespace QuadP1P1Cell
  }  // namespace UFEM
}  // namespace COOLFluiD

