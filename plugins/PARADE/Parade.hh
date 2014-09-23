#ifndef COOLFluiD_Parade_hh
#define COOLFluiD_Parade_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement computation of aerodynamic coefficients
  namespace Parade {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Parade
 */
class ParadeModule : public Environment::ModuleRegister<ParadeModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Parade";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an interface for the PARADE library.";
  }

}; // end ParadeModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace Parade

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Parade_hh

