#ifndef COOLFluiD_SpectralFV_hh
#define COOLFluiD_SpectralFV_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement Spectral Finite Volume space discretization
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module SpectralFV
 */
class SpectralFVModule : public Environment::ModuleRegister<SpectralFVModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "SpectralFV";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the Spectral Finite Volume space discretization.";
  }

}; // end SpectralFVModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_SpectralFV_hh

