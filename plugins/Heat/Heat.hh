#ifndef COOLFluiD_Physics_Heat_hh
#define COOLFluiD_Physics_Heat_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement Heat transfer physical model.
  namespace Heat {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Heat
 */
class HeatModule : public Environment::ModuleRegister<HeatModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Heat";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a Heat transfer physical model.";
  }

}; // end HeatModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_Heat_hh
