#ifndef COOLFluiD_Physics_Poisson_hh
#define COOLFluiD_Physics_Poisson_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a Poisson equation physical model.
  namespace Poisson {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Poisson
 */
class PoissonModule : public Environment::ModuleRegister<PoissonModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Poisson";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a Poisson equation physical model.";
  }

}; // end PoissonModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_Poisson_hh
