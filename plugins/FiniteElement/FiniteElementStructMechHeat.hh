#ifndef COOLFluiD_Numerics_FiniteElementStructMechHeat_hh
#define COOLFluiD_Numerics_FiniteElementStructMechHeat_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement the FiniteElementStructMechHeat space discretization method
  namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteElementStructMechHeat
 */
class FiniteElementStructMechHeatModule : public Environment::ModuleRegister<FiniteElementStructMechHeatModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteElementStructMechHeat";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the StructMechHeat bindings for the FiniteElement space discretization method.";
  }

}; // end FiniteElementStructMechHeatModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteElementStructMechHeat_hh

