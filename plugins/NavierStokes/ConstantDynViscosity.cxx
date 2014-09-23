#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"

#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes/ConstantDynViscosity.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ConstantDynViscosity,
                            DynamicViscosityLaw,
                            NavierStokesModule,
                            1>
constantDynViscosityProvider("Constant");
      
//////////////////////////////////////////////////////////////////////////////

void ConstantDynViscosity::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("Value", "Constant value");
}

//////////////////////////////////////////////////////////////////////////////

ConstantDynViscosity::ConstantDynViscosity(const std::string& name) :
DynamicViscosityLaw(name)
{
  addConfigOptionsTo(this);
  
  m_value = 0.0;
  setParameter("Value", &m_value);
}
      
//////////////////////////////////////////////////////////////////////////////

ConstantDynViscosity::~ConstantDynViscosity()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConstantDynViscosity::configure ( Config::ConfigArgs& args )
{
  DynamicViscosityLaw::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConstantDynViscosity::compute(const CFreal& pdim, const CFreal& Tdim)
{
  return m_value;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
