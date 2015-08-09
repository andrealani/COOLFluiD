#ifndef COOLFluiD_Physics_Plato_hh
#define COOLFluiD_Physics_Plato_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    /// The classes that implement an interface to the Plato library.
    namespace Plato {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Plato
 */
class PlatoModule : public Environment::ModuleRegister<PlatoModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Plato";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an interface to the Plato library.";
  }

}; // end PlatoModule

//////////////////////////////////////////////////////////////////////////////

    }  // namespace Plato

  }  // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Plato_hh
