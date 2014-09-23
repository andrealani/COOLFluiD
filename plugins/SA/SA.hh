#ifndef COOLFluiD_Physics_SA_hh
#define COOLFluiD_Physics_SA_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement Spalart-Almaras flow physical model.
  namespace SA {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module SA
 */
class SAModule : public Environment::ModuleRegister<SAModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "SA";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a Spalart-Almaras flow physical model.";
  }

}; // end SAModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_SA_hh
