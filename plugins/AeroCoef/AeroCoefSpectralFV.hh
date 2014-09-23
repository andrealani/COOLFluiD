#ifndef COOLFluiD_Numerics_AeroCoefSpectralFV_hh
#define COOLFluiD_Numerics_AeroCoefSpectralFV_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement computation of aerodynamic coefficients
    /// with bindings to spectral volume methods
  namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module AeroCoefSpectralFV
 */
class AeroCoefSpectralFVModule : public Environment::ModuleRegister<AeroCoefSpectralFVModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "AeroCoefSpectralFV";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements computation of aerodynamic coefficients with bindings to spectral difference methods.";
  }

}; // end AeroCoefSpectralFVModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_AeroCoefSpectralFV_hh
