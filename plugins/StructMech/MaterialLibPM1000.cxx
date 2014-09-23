
#include "Environment/ObjectProvider.hh"

#include "StructMech/StructMech.hh"
#include "StructMech/MaterialLibPM1000.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaterialLibPM1000,
                            MaterialPropertyLib,
                            StructMechModule,
                            1>
MaterialLibPM1000Provider("PM1000");

//////////////////////////////////////////////////////////////////////////////

void MaterialLibPM1000::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("YoungModulus", "Young Modulus");
  options.addConfigOption< CFreal > ("Poisson", "Poisson Coeficient");
  options.addConfigOption< CFreal > ("Density", "Density");
  options.addConfigOption< CFreal > ("ThermalExpansion", "Thermal Expansion Coef");
}

//////////////////////////////////////////////////////////////////////////////

MaterialLibPM1000::MaterialLibPM1000(const std::string& name) :
MaterialPropertyLib(name)
{
  addConfigOptionsTo(this);

  //these 2 should vary with temperature (Values at 1273 kelvin)
  m_alpha     = 0.000017;
  m_young     = 119E9;

  m_poisson   = 0.30;
  m_density   = 8240.;
  m_isAnisotropic = false;

  setParameter("YoungModulus", &m_young);
  setParameter("Poisson", &m_poisson);
  setParameter("Density", &m_density);
  setParameter("ThermalExpansion", &m_alpha);

}

//////////////////////////////////////////////////////////////////////////////

MaterialLibPM1000::~MaterialLibPM1000()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibPM1000::computeYoungModulus()
{
  return m_young;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibPM1000::computePoissonCoef()
{
  return m_poisson;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibPM1000::computeThermalExpansionCoef()
{
  return m_alpha;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
