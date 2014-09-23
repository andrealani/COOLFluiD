#ifndef COOLFluiD_Physics_Mutation2_hh
#define COOLFluiD_Physics_Mutation2_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    /// The classes that implement an interface to the Mutation2 library.
    namespace Mutation2 {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Mutation2
 */
class Mutation2Module : public Environment::ModuleRegister<Mutation2Module> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Mutation2";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an interface to the Mutation2 library.";
  }

}; // end Mutation2Module

//////////////////////////////////////////////////////////////////////////////

    }  // namespace Mutation2

  }  // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Mutation2_hh
