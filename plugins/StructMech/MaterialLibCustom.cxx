
#include "Environment/ObjectProvider.hh"

#include "StructMech/StructMech.hh"
#include "StructMech/MaterialLibCustom.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaterialLibCustom,
                            MaterialPropertyLib,
                            StructMechModule,
                            1>
MaterialLibCustomProvider("Custom");

//////////////////////////////////////////////////////////////////////////////

void MaterialLibCustom::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("YoungModulus", "Young Modulus");
  options.addConfigOption< CFreal > ("Poisson", "Poisson Coeficient");
  options.addConfigOption< CFreal > ("Density", "Density");
  options.addConfigOption< CFreal > ("ThermalExpansion", "Thermal Expansion Coef");
}

//////////////////////////////////////////////////////////////////////////////

MaterialLibCustom::MaterialLibCustom(const std::string& name) :
MaterialPropertyLib(name)
{
  addConfigOptionsTo(this);

  // defaults are from
  m_young     = 200E9;
  m_poisson   = 0.33;
  m_density   = 7850.;
  m_alpha     = 0.000012;
  m_isAnisotropic = false;

  setParameter("YoungModulus", &m_young);
  setParameter("Poisson", &m_poisson);
  setParameter("Density", &m_density);
  setParameter("ThermalExpansion", &m_alpha);

}

//////////////////////////////////////////////////////////////////////////////

MaterialLibCustom::~MaterialLibCustom()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibCustom::computeYoungModulus()
{
  return m_young;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibCustom::computePoissonCoef()
{
  return m_poisson;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibCustom::computeThermalExpansionCoef()
{
  return m_alpha;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
