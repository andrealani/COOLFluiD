#ifndef COOLFluiD_Numerics_FiniteElementHeat_hh
#define COOLFluiD_Numerics_FiniteElementHeat_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement the FiniteElementHeat space discretization method
  namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteElementHeat
 */
class FiniteElementHeatModule : public Environment::ModuleRegister<FiniteElementHeatModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteElementHeat";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the Heat bindings for the FiniteElement space discretization method.";
  }

}; // end FiniteElementHeatModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteElementHeat_hh

