#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"

#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes/SimplerSutherlandDynViscosity.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SimplerSutherlandDynViscosity,
                            DynamicViscosityLaw,
                            NavierStokesModule,
                            1>
SimplerSutherlandDynViscosityProvider("SimplerSutherland");

//////////////////////////////////////////////////////////////////////////////

void SimplerSutherlandDynViscosity::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string > ("Gas", "Name of the gas to model");
  options.addConfigOption< CFreal >   ("ViscRef", "Reference viscosity");
  options.addConfigOption< CFreal >   ("TRef", "Reference temperature");
}

//////////////////////////////////////////////////////////////////////////////

SimplerSutherlandDynViscosity::SimplerSutherlandDynViscosity(const std::string& name) :
DynamicViscosityLaw(name)
{
  addConfigOptionsTo(this);
  
  m_gas = "";

  // defaults are for air
  m_TRef      = 110.4;
  m_ViscRef   = 0.000001458;

  setParameter("TRef",    &m_TRef);
  setParameter("ViscRef", &m_ViscRef);
  setParameter("Gas",     &m_gas);
}

//////////////////////////////////////////////////////////////////////////////

SimplerSutherlandDynViscosity::~SimplerSutherlandDynViscosity()
{
}

//////////////////////////////////////////////////////////////////////////////

void SimplerSutherlandDynViscosity::configure ( Config::ConfigArgs& args )
{
  DynamicViscosityLaw::configure(args);

  if ( !m_gas.empty() )
  {
    if ( m_gas == "Air")
    {
      m_TRef      = 110.4;
      m_ViscRef   = 0.000001458;
    }
    if ( m_gas == "CO")
    {
      m_TRef      = 109.;
      m_ViscRef   = 1.40E-6;
    }
    if ( m_gas == "CO2")
    {
      m_TRef      = 233.;
      m_ViscRef   = 1.55E-6;
    }
    if ( m_gas == "N2")
    {
      m_TRef      = 102.;
      m_ViscRef   = 1.39E-6;
    }
    if ( m_gas == "H2")
    {
      m_TRef      = 71.;
      m_ViscRef   = 0.65E-6;
    }
    if ( m_gas == "O2")
    {
      m_TRef      = 110.;
      m_ViscRef   = 1.65E-6;
    }
    if ( ! ( m_gas == "Air" ||  m_gas == "CO" || m_gas == "CO2" || m_gas == "N2"  || m_gas == "H2"  || m_gas == "O2" ) )
    {
      throw Common::BadValueException (FromHere(),"Viscosity unknown for gas : " + m_gas);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal SimplerSutherlandDynViscosity::compute(const CFreal& pdim, const CFreal& Tdim)
{
  return m_ViscRef*std::pow(Tdim,1.5)/( m_TRef + Tdim );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
