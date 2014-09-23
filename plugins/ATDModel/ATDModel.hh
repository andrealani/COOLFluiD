#ifndef COOLFluiD_Physics_ATDModel_hh
#define COOLFluiD_Physics_ATDModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    /// The classes that implement an interface to the ATDModel library.
    namespace ATDModel {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module ATDModel
 */
class ATDModelModule : public Environment::ModuleRegister<ATDModelModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "ATDModel";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an interface to the ATDModel library.";
  }

}; // end ATDModelModule

//////////////////////////////////////////////////////////////////////////////

    }  // namespace ATDModel

  }  // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ATDModel_hh
