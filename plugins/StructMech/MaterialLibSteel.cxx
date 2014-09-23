
#include "Environment/ObjectProvider.hh"

#include "StructMech/StructMech.hh"
#include "StructMech/MaterialLibSteel.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaterialLibSteel,
                            MaterialPropertyLib,
                            StructMechModule,
                            1>
MaterialLibSteelProvider("Steel");

//////////////////////////////////////////////////////////////////////////////

void MaterialLibSteel::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("YoungModulus", "Young Modulus");
  options.addConfigOption< CFreal > ("Poisson", "Poisson Coeficient");
  options.addConfigOption< CFreal > ("Density", "Density");
  options.addConfigOption< CFreal > ("ThermalExpansion", "Thermal Expansion Coef");
}

//////////////////////////////////////////////////////////////////////////////

MaterialLibSteel::MaterialLibSteel(const std::string& name) :
MaterialPropertyLib(name)
{
  addConfigOptionsTo(this);

  // defaults are from
  m_young     = 205e9;
  m_poisson   = 0.3;
  m_density   = 7850.;
  m_alpha     = 0.000012;
  m_isAnisotropic = false;

  setParameter("YoungModulus", &m_young);
  setParameter("Poisson", &m_poisson);
  setParameter("Density", &m_density);
  setParameter("ThermalExpansion", &m_alpha);
}

//////////////////////////////////////////////////////////////////////////////

MaterialLibSteel::~MaterialLibSteel()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibSteel::computeYoungModulus()
{
  return m_young;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibSteel::computePoissonCoef()
{
  return m_poisson;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibSteel::computeThermalExpansionCoef()
{
  return m_alpha;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibSteel::computeDensity()
{
  return m_density;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
