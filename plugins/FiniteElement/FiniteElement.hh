#ifndef COOLFluiD_Numerics_FiniteElement_hh
#define COOLFluiD_Numerics_FiniteElement_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement the FiniteElement space discretization method
  namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteElement
 */
class FiniteElementModule : public Environment::ModuleRegister<FiniteElementModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteElement";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the FiniteElement space discretization method.";
  }

}; // end FiniteElementModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteElement_hh

