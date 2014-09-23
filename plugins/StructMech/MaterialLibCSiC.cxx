
#include "Environment/ObjectProvider.hh"

#include "StructMech/StructMech.hh"
#include "StructMech/MaterialLibCSiC.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaterialLibCSiC,
                            MaterialPropertyLib,
                            StructMechModule,
                            1>
MaterialLibCSiCProvider("CSiC");

//////////////////////////////////////////////////////////////////////////////

void MaterialLibCSiC::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("YoungModulus", "Young Modulus");
  options.addConfigOption< CFreal > ("Poisson", "Poisson Coeficient");
  options.addConfigOption< CFreal > ("Density", "Density");
  options.addConfigOption< CFreal > ("ThermalExpansion", "Thermal Expansion Coef");
}

//////////////////////////////////////////////////////////////////////////////

MaterialLibCSiC::MaterialLibCSiC(const std::string& name) :
MaterialPropertyLib(name)
{
  addConfigOptionsTo(this);

  // defaults are from http://www.accuratus.com/silicar.html
  m_young     = 410E9;
  m_poisson   = 0.14;
  m_density   = 3100.;
  m_alpha     = 0.000004;
  m_isAnisotropic = false;

  setParameter("YoungModulus", &m_young);
  setParameter("Poisson", &m_poisson);
  setParameter("Density", &m_density);
  setParameter("ThermalExpansion", &m_alpha);

}

//////////////////////////////////////////////////////////////////////////////

MaterialLibCSiC::~MaterialLibCSiC()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibCSiC::computeYoungModulus()
{
  return m_young;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibCSiC::computePoissonCoef()
{
  return m_poisson;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibCSiC::computeThermalExpansionCoef()
{
  return m_alpha;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
