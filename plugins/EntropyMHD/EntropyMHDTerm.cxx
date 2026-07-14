#include "EntropyMHD/EntropyMHDTerm.hh"
#include "Config/ConfigObject.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

void EntropyMHDTerm::defineConfigOptions(Config::OptionList& options)
{
  // gamma
  options.addConfigOption< CFreal >("gamma", "Ratio of specific heats (5/3 for thermodynamic).");

  // projection scheme
  options.addConfigOption< CFreal, Config::DynamicOption<> >("refSpeed", "Reference speed for projection scheme.");

  // reference values for non-dimensionalization
  options.addConfigOption< CFreal >("nRef", "Reference proton density [m^-3].");
  options.addConfigOption< CFreal >("BRef", "Reference magnetic field [T].");
  options.addConfigOption< CFreal >("lRef", "Reference length (R_sun) [m].");
  options.addConfigOption< CFreal >("TRef", "Reference temperature [K].");
  options.addConfigOption< CFreal >("mass", "Mass of central object [kg].");

  // Thermodynamic parameters
  options.addConfigOption< CFreal >("SpitzerKappa0", "Spitzer conductivity kappa_0 [W m^-1 K^-7/2].");
  options.addConfigOption< CFreal >("HeatingH0", "Coronal heating amplitude H_0 [SI].");
  options.addConfigOption< CFreal >("HeatingLambda", "Coronal heating scale height [R_sun].");
  options.addConfigOption< CFreal >("TRBroadeningTcut", "Transition region broadening cutoff temperature [K].");

}

//////////////////////////////////////////////////////////////////////////////

EntropyMHDTerm::EntropyMHDTerm(const std::string& name)
  : BaseTerm(name)
{
  addConfigOptionsTo(this);

  // default gamma = 5/3 (thermodynamic)
  _gamma = 5.0/3.0;
  setParameter("gamma", &_gamma);

  _refSpeed = 1.0;
  setParameter("refSpeed", &_refSpeed);

  // Reference values (COCONUT defaults)
  _nRef = 1.0e14;
  setParameter("nRef", &_nRef);

  _BRef = 2.2e-4;
  setParameter("BRef", &_BRef);

  _lRef = 6.9551e8;
  setParameter("lRef", &_lRef);

  _TRef = 1.5e6;
  setParameter("TRef", &_TRef);

  _mass = 1.989e30;
  setParameter("mass", &_mass);

  // Thermodynamic defaults
  _kappa0 = 9.0e-12;
  setParameter("SpitzerKappa0", &_kappa0);

  _heatingH0 = 4.0e-2;
  setParameter("HeatingH0", &_heatingH0);

  _heatingLambda = 0.7;
  setParameter("HeatingLambda", &_heatingLambda);

  _TRBroadeningTcut = 5.0e5;
  setParameter("TRBroadeningTcut", &_TRBroadeningTcut);

  _gravFactorNonDim = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

EntropyMHDTerm::~EntropyMHDTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHDTerm::configure(Config::ConfigArgs& args)
{
  BaseTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHDTerm::setupPhysicalData()
{
  cf_assert(getDataSize() == 14);

  // Compute non-dimensional gravity factor: g_dimless = -gravFactor / r^2
  const CFreal GMsun = 1.327474512e20;  // m^3 s^-2
  const CFreal mu0   = 1.2566370614e-6; // H/m
  const CFreal mp    = 1.67262158e-27;
  const CFreal me    = 9.10938188e-31;
  const CFreal rhoRef = _nRef * (mp + me);
  const CFreal V0 = _BRef / std::sqrt(mu0 * rhoRef);
  const CFreal g0 = V0 * V0 / _lRef;
  _gravFactorNonDim = GMsun / (_lRef * _lRef * g0);

  CFLog(INFO, "EntropyMHDTerm: gravFactor (non-dim) = " << _gravFactorNonDim
        << ", V0 = " << V0 << " m/s\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
