#ifndef COOLFluiD_Physics_Mutation_hh
#define COOLFluiD_Physics_Mutation_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    /// The classes that implement an interface to the Mutation library.
    namespace Mutation {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Mutation
 */
class MutationModule : public Environment::ModuleRegister<MutationModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Mutation";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an interface to the Mutation library.";
  }

}; // end MutationModule

//////////////////////////////////////////////////////////////////////////////

    }  // namespace Mutation

  }  // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Mutation_hh
