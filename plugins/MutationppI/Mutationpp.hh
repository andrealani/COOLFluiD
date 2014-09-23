#ifndef COOLFluiD_Physics_Mutationpp_hh
#define COOLFluiD_Physics_Mutationpp_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    /// The classes that implement an interface to the Mutationpp library.
    namespace Mutationpp {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Mutationpp
 */
class MutationppModule : public Environment::ModuleRegister<MutationppModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Mutationpp";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an interface to the Mutationpp library.";
  }

}; // end MutationppModule

//////////////////////////////////////////////////////////////////////////////

    }  // namespace Mutationpp

  }  // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Mutationpp_hh
