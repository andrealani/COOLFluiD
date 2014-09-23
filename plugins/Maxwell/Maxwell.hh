#ifndef COOLFluiD_Physics_Maxwell_hh
#define COOLFluiD_Physics_Maxwell_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement Maxwell physical model.
  namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Maxwell
 */
class MaxwellModule : public Environment::ModuleRegister<MaxwellModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Maxwell";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a Maxwell physical model.";
  }

}; // end MaxwellModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_Maxwell_hh
