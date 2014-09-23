#ifndef COOLFluiD_Numerics_FiniteElementStructMech_hh
#define COOLFluiD_Numerics_FiniteElementStructMech_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement the FiniteElementStructMech space discretization method
  namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteElementStructMech
 */
class FiniteElementStructMechModule : public Environment::ModuleRegister<FiniteElementStructMechModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteElementStructMech";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the StructMech bindings for the FiniteElement space discretization method.";
  }

}; // end FiniteElementStructMechModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteElementStructMech_hh

