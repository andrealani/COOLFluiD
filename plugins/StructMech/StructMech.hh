#ifndef COOLFluiD_Physics_StructMech_hh
#define COOLFluiD_Physics_StructMech_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a structural mechanics physical model.
  namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module StructMech
 */
class StructMechModule : public Environment::ModuleRegister<StructMechModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "StructMech";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a structural mechanics physical model.";
  }

}; // end StructMechModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_StructMech_hh
