#ifndef COOLFluiD_Physics_KOmega_hh
#define COOLFluiD_Physics_KOmega_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a KOmega flow physical model.
  namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module KOmega
 */
class KOmegaModule : public Environment::ModuleRegister<KOmegaModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "KOmega";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a KOmega flow physical model.";
  }

}; // end KOmegaModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_KOmega_hh
