#ifndef COOLFluiD_Physics_NEQ_hh
#define COOLFluiD_Physics_NEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a NEQ flow physical model.
  namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module NEQ
 */
class NEQModule : public Environment::ModuleRegister<NEQModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "NEQ";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a NEQ flow physical model.";
  }

}; // end NEQModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_NEQ_hh
