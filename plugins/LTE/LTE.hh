#ifndef COOLFluiD_Physics_LTE_hh
#define COOLFluiD_Physics_LTE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a LTE flow physical model.
  namespace LTE {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module LTE
 */
class LTEModule : public Environment::ModuleRegister<LTEModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "LTE";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a LTE flow physical model.";
  }

}; // end LTEModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_LTE_hh
