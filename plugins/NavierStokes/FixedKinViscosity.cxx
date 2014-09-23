

#include "Environment/ObjectProvider.hh"

#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes/FixedKinViscosity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FixedKinViscosity,
                            DynamicViscosityLaw,
                            NavierStokesModule,
                            1>
fixedKinematicViscosityProvider("FixedKinematicViscosity");

//////////////////////////////////////////////////////////////////////////////

void FixedKinViscosity::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("KinVisc", "Value of fixed kinematic viscosity");
}

//////////////////////////////////////////////////////////////////////////////

FixedKinViscosity::FixedKinViscosity(const std::string& name) :
DynamicViscosityLaw(name)
{
  addConfigOptionsTo(this);

  m_fixedKinVisc = 0.0;
  setParameter("KinVisc", &m_fixedKinVisc);
}

//////////////////////////////////////////////////////////////////////////////

FixedKinViscosity::~FixedKinViscosity()
{
}
//////////////////////////////////////////////////////////////////////////////

CFreal FixedKinViscosity::compute(const CFreal& pdim, const CFreal& Tdim)
{
//  return pdim/(m_eulerVarSet->getModel()->getRdim()*Tdim)*m_fixedKinVisc;
  return pdim/(287.046*Tdim)*m_fixedKinVisc;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
