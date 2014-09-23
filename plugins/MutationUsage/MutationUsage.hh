#ifndef COOLFluiD_MutationUsage_hh
#define COOLFluiD_MutationUsage_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement an application to use Mutation library for computing physical coefficients
  namespace MutationUsage {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module MutationUsage
 */
class MutationUsageModule : public Environment::ModuleRegister<MutationUsageModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "MutationUsage";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an application to use Mutation library for computing physical coefficients.";
  }

}; // end MutationUsageModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace MutationUsage

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_MutationUsage_hh
