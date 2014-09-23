#ifndef COOLFluiD_Physics_Chemistry_CH4_hh
#define COOLFluiD_Physics_Chemistry_CH4_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

  namespace Chemistry {

    /// The classes that implement CH4 physical model.
  namespace CH4 {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module CH4
 */
class CH4Module : public Environment::ModuleRegister<CH4Module> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "CH4";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a Chemistry CH4 equation Physical model";
  }

}; // end CH4Module

//////////////////////////////////////////////////////////////////////////////

    } // namespace CH4

    }  // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Chemistry_CH4_hh

