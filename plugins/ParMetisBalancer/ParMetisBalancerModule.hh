#ifndef COOLFluiD_Numerics_ParMetisBalancerModule_hh
#define COOLFluiD_Numerics_ParMetisBalancerModule_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    /// The classes that implement Dynamic Balancing
    namespace ParMetisBalancer {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module 
 */
class ParMetisBalancerModule : public Environment::ModuleRegister<ParMetisBalancerModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "ParMetisBalancer";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implementes the Dynamic Balancing through use of ParMetis.";
  }

}; // end ParMetisBalancerModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace ParMetisBalancer



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_hh
