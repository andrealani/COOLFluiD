#ifndef COOLFluiD_Physics_StructMechHeat_StructMechHeat_hh
#define COOLFluiD_Physics_StructMechHeat_StructMechHeat_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a structural mechanics physical model.
  namespace StructMechHeat {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module StructMechHeat
 */
class StructMechHeatModule : public Environment::ModuleRegister<StructMechHeatModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "StructMechHeat";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a structural mechanics physical model including temperature effect.";
  }

}; // end StructMechHeatModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_StructMechHeat_StructMechHeat_hh
