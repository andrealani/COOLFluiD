
#include "Environment/ObjectProvider.hh"

#include "Common/BadValueException.hh"

#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes/SutherlandDynViscosity.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SutherlandDynViscosity,
                            DynamicViscosityLaw,
                            NavierStokesModule,
                            1>
SutherlandDynViscosityProvider("Sutherland");

//////////////////////////////////////////////////////////////////////////////

void SutherlandDynViscosity::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string > ("Gas", "Name of the gas to model");
  options.addConfigOption< CFreal >   ("SuthConst", "Sutherland constant for the fluid");
  options.addConfigOption< CFreal >   ("ViscRef", "Reference viscosity");
  options.addConfigOption< CFreal >   ("TRef", "Reference temperature at which the reference viscosity is measured");
}

//////////////////////////////////////////////////////////////////////////////

SutherlandDynViscosity::SutherlandDynViscosity(const std::string& name) :
DynamicViscosityLaw(name)
{
  addConfigOptionsTo(this);

  // defaults are for air
  m_SuthConst = 120.;
  m_TRef      = 291.15;
  m_ViscRef   = 18.27E-6;

  setParameter("SuthConst", &m_SuthConst);
  setParameter("ViscRef",   &m_ViscRef);
  setParameter("TRef",      &m_TRef);
  setParameter("Gas",       &m_gas);
}

//////////////////////////////////////////////////////////////////////////////

SutherlandDynViscosity::~SutherlandDynViscosity()
{
}

//////////////////////////////////////////////////////////////////////////////

void SutherlandDynViscosity::configure ( Config::ConfigArgs& args )
{
  DynamicViscosityLaw::configure(args);

  if ( !m_gas.empty() )
  {
    if ( m_gas == "Air")
    {
      m_SuthConst = 120.;
      m_TRef      = 291.15;
      m_ViscRef   = 18.27E-6;
    }
    if ( m_gas == "CO")
    {
      m_SuthConst = 118.;
      m_TRef      = 288.15;
      m_ViscRef   = 17.2E-6;
    }
    if ( m_gas == "CO2")
    {
      m_SuthConst = 240.;
      m_TRef      = 293.15;
      m_ViscRef   = 14.8E-6;
    }
    if ( m_gas == "H2")
    {
      m_SuthConst = 72.;
      m_TRef      = 293.85;
      m_ViscRef   = 8.76E-6;
    }
    if ( m_gas == "O2")
    {
      m_SuthConst = 127.;
      m_TRef      = 292.25;
      m_ViscRef   = 20.18E-6;
    }
    if ( ! ( m_gas == "Air" ||  m_gas == "CO" || m_gas == "CO2" || m_gas == "H2" || m_gas == "O2" ) )
    {
      throw Common::BadValueException (FromHere(),"Viscosity unknown for gas : " + m_gas);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal SutherlandDynViscosity::compute(const CFreal& pdim, const CFreal& Tdim)
{
  return m_ViscRef*(m_TRef+m_SuthConst)/(Tdim+m_SuthConst)*std::pow(Tdim/m_TRef,1.5);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
