

#include "Environment/ObjectProvider.hh"

#include "StructMech/StructMech.hh"
#include "StructMech/MaterialLibHTI.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaterialLibHTI,
                            MaterialPropertyLib,
                            StructMechModule,
                            1>
MaterialLibHTIProvider("HTI");

//////////////////////////////////////////////////////////////////////////////

void MaterialLibHTI::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("YoungModulus", "Young Modulus");
  options.addConfigOption< CFreal > ("Poisson", "Poisson Coeficient");
  options.addConfigOption< CFreal > ("Density", "Density");
  options.addConfigOption< CFreal > ("ThermalExpansion", "Thermal Expansion Coef");
}

//////////////////////////////////////////////////////////////////////////////

MaterialLibHTI::MaterialLibHTI(const std::string& name) :
MaterialPropertyLib(name)
{
  addConfigOptionsTo(this);

  // defaults are from N/A (small values)
  m_young     = 10000000E9;
  m_poisson   = 0.;
  m_density   = 1.;
  m_alpha     = 0.0;
  m_isAnisotropic = false;

  setParameter("YoungModulus", &m_young);
  setParameter("Poisson", &m_poisson);
  setParameter("Density", &m_density);
  setParameter("ThermalExpansion", &m_alpha);

}

//////////////////////////////////////////////////////////////////////////////

MaterialLibHTI::~MaterialLibHTI()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibHTI::computeYoungModulus()
{
  return m_young;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibHTI::computePoissonCoef()
{
  return m_poisson;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibHTI::computeThermalExpansionCoef()
{
  return m_alpha;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
