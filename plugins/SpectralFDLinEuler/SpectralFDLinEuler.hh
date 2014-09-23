#ifndef COOLFluiD_SpectralFDLinEuler_hh
#define COOLFluiD_SpectralFDLinEuler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the Spectral Finite Difference space discretization.
  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module SpectralFDLinEuler
 */
class SpectralFDLinEulerModule : public Environment::ModuleRegister<SpectralFDLinEulerModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "SpectralFDLinEuler";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the classes specific for the Spectral Finite Difference space discretization for linearized Euler.";
  }

}; // end SpectralFDLinEulerModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_SpectralFDLinEuler_hh

