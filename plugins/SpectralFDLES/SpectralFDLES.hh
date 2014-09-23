#ifndef COOLFluiD_SpectralFDLES_hh
#define COOLFluiD_SpectralFDLES_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the Spectral Finite Difference space discretization.
  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module SpectralFDLES
 * @author Kris Van den Abeele
 * @author Matteo Parsani
 */
class SpectralFDLESModule : public Environment::ModuleRegister<SpectralFDLESModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "SpectralFDLES";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the classes specific for the Spectral Finite Difference space discretization for LES.";
  }

}; // end SpectralFDLESModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_SpectralFDLES_hh

