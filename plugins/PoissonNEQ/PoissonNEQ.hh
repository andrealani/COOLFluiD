#ifndef COOLFluiD_Physics_PoissonNEQ_hh
#define COOLFluiD_Physics_PoissonNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

  namespace PoissonNEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module PoissonNEQ
 */
class PoissonNEQModule : public Environment::ModuleRegister<PoissonNEQModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "PoissonNEQ";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a PoissonNEQ flow physical model.";
  }

}; // end PoissonNEQModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace PoissonNEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_PoissonNEQ_hh
