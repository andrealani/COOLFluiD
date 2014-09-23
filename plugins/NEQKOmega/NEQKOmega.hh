#ifndef COOLFluiD_Physics_NEQKOmega_hh
#define COOLFluiD_Physics_NEQKOmega_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a NEQKOmega flow physical model.
  namespace NEQKOmega {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module NEQKOmega
 */
class NEQKOmegaModule : public Environment::ModuleRegister<NEQKOmegaModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "NEQKOmega";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a NEQKOmega flow physical model.";
  }

}; // end NEQKOmegaModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQKOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_NEQKOmega_hh
