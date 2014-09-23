#ifndef COOLFluiD_Numerics_FiniteElementChemistry_hh
#define COOLFluiD_Numerics_FiniteElementChemistry_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement the FiniteElementChemistry space discretization method
  namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteElementChemistry
 */
class FiniteElementChemistryModule : public Environment::ModuleRegister<FiniteElementChemistryModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteElementChemistry";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the Chemistry bindings for the FiniteElement space discretization method.";
  }

}; // end FiniteElementChemistryModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteElementChemistry_hh

