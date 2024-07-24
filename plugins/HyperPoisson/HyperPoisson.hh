#ifndef COOLFluiD_Physics_HyperPoisson_hh
#define COOLFluiD_Physics_HyperPoisson_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace HyperPoisson {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module HyperPoisson
 */
class HyperPoissonModule : public Environment::ModuleRegister<HyperPoissonModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "HyperPoisson";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a Hyperbolized Poisson physical model.";
  }

}; // end HyperPoissonModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace HyperPoisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_HyperPoisson_hh
