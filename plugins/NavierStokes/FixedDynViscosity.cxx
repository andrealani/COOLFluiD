

#include "Environment/ObjectProvider.hh"

#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes/FixedDynViscosity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FixedDynViscosity,
                            DynamicViscosityLaw,
                            NavierStokesModule,
                            1>
fixedDynViscosityProvider("Fixed");

//////////////////////////////////////////////////////////////////////////////

void FixedDynViscosity::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("Fix_visc", "Fixed viscosity choosed");
}

//////////////////////////////////////////////////////////////////////////////

FixedDynViscosity::FixedDynViscosity(const std::string& name) :
DynamicViscosityLaw(name)
{
  addConfigOptionsTo(this);

  // defaults are the one of Sutherland for P and T atmospheric
   m_FixedVisc = 1.7e-5;

  setParameter("Fix_visc", &m_FixedVisc);

}

//////////////////////////////////////////////////////////////////////////////

FixedDynViscosity::~FixedDynViscosity()
{
}
//////////////////////////////////////////////////////////////////////////////

CFreal FixedDynViscosity::compute(const CFreal& pdim, const CFreal& Tdim)
{
  return m_FixedVisc;
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
