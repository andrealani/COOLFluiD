#ifndef COOLFluiD_SpectralFD_hh
#define COOLFluiD_SpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement Spectral Finite Difference space discretization
  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module SpectralFD
 */
class SpectralFDModule : public Environment::ModuleRegister<SpectralFDModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "SpectralFD";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the Spectral Finite Difference space discretization.";
  }

}; // end SpectralFDModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_SpectralFD_hh

