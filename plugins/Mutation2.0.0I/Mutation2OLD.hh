#ifndef COOLFluiD_Physics_Mutation2OLD_hh
#define COOLFluiD_Physics_Mutation2OLD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    /// The classes that implement an interface to the Mutation2 library.
    namespace Mutation2OLD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Mutation2OLD
 */
class Mutation2OLDModule : public Environment::ModuleRegister<Mutation2OLDModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Mutation2OLD";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an interface to the Mutation2OLD library.";
  }

}; // end Mutation2OLDModule

//////////////////////////////////////////////////////////////////////////////

    }  // namespace Mutation2OLD

  }  // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Mutation2OLD_hh
